import numpy as np
import h5py
import scipy
import json
import os

from ossdbs.electrodes.defaults import default_electrode_parameters
from ossdbs.utils.settings import Settings

class LeadSettings:
    
    def __init__(self, mat_file_path):
        try:
            self._file = h5py.File(str(mat_file_path), 'r')
            self._h5 = True
        except:
            self._file = scipy.io.loadmat(mat_file_path)
            self._h5 = False
        self._settings = self._file["settings"]
        
        # Check to see that both electrodes are either current controlled or voltage controlled
        #if file.root.settings.current_control[0][0]!=file.root.settings.current_control[0][1]:
        cur_ctrl_arr = self.get_cur_ctrl()
        if all(~np.isnan(cur_ctrl_arr)):
            if cur_ctrl_arr[0] != cur_ctrl_arr[1]:
                print("Simultaneous use of VC and CC is not allowed for safety reasons!")
                raise SystemExit
        
        self.NUM_ELECS = 1
    
    
    # Make the lead settings into an OSS-DBS json 
    def make_oss_settings(self, hemis_idx=0):

        #index_side = hemis_idx
        hemis_idx = int(hemis_idx)

        # Constants
        # OSS strings used for tissue types 
        UNK = "Unknown"
        CSF = "CSF"
        GM = "Gray matter"
        WM = "White matter"
        BLOOD = "Blood"
        # LEAD to OSS tissue name map
        TISSUE_NAMES = {"Unknown":UNK,"CSF":CSF,"GM":GM,"WM":WM,"Blood":BLOOD}
        HEMIS_OUTPUT_PATHS = ["Results_rh","Results_lh"]
        
        oss_elec_name = self.get_elec_type().replace(" ", "")
        # Check that oss electrode name is valid
        if oss_elec_name not in default_electrode_parameters.keys():
            raise Exception(oss_elec_name + " is not a recognized electrode type")
        
        
        ### CALCULATE TIP POSITION FROM IMPLANTATION COORDINATE AND SECOND COORDINATE ###
        
        # Find offset between lead implantation coordinate at center of contact
        # and the tip position of the electrode
        elec_params = default_electrode_parameters[oss_elec_name]
        imp_coords = np.array(self.get_imp_coord())
        sec_coords = np.array(self.get_sec_coord())
        directions = sec_coords - imp_coords
        unit_directions = directions/np.linalg.norm(directions,axis=1)[:,np.newaxis]
        tip_pos = imp_coords - elec_params.offset * unit_directions
        
        
        ### SORT OUT CURRENT V VOLTAGE DETAILS ###
        
        pulse_amps = self.get_phi_vec()

        # Set either voltage or current to pulse_amp values based on
        # whether it's current controlled or voltage controlled
        pulse_amp_key = "Voltage[V]"
        not_pulse_amp_key = "Current[A]"
        # Assuming electrodes are both current controlled or both voltage controlled
        if self.get_cur_ctrl()[hemis_idx]:
            temp = pulse_amp_key
            pulse_amp_key = not_pulse_amp_key
            not_pulse_amp_key = temp
            # Lead-DBS uses mA as the input
            pulse_amps = pulse_amps * 0.001
        # Fix of VC random grounding bug for Lead-DBS stim settings
        else:
            pulse_amps[pulse_amps == 0] = float('nan')
        
        # mono_elec: whether or not one of the electrodes in the simulation is monopolar
        # This is used for determining whether there needs to be a non-electrode ground
        mono_elec = False
        
        
        ### MAKE THE ITERABLE ELECTRODE AND CONTACT SUB-DICTIONARIES ###
        
        # Make the dictionaries for the electrode settings one by one
        elec_dicts = np.empty(self.get_num_elecs(), dtype=object)
        for index_side in [hemis_idx]:
            
            pulse_amp = pulse_amps[index_side,:]
            
            # Determining the value of mono_elec for voltage_controlled
            # If you haven't found a monopolar electrode yet, keep looking
            # Once you have found one, you don't need to keep checking to
            # see if there is one
            if not mono_elec:
                # A monopolar elec in cur ctrl has a pulse_amp array that doesn't sum to 0
                if self.get_cur_ctrl()[index_side]:
                    mono_elec = np.sum(pulse_amp[~np.isnan(pulse_amp)]) != 0
                # A monopolar elec in volt ctrl has a pulse_amp with two unique values: nan and x
                else:
                    mono_elec = bool(self.get_case_grnd()[index_side])
                
            # cntct_dicts is a list of the contacts that will go in the json for this electrode
            # We only include the active ones
            cntct_dicts = np.empty(np.sum(~np.isnan(pulse_amp)),dtype=object)
            cntcts_made = 0

            for i in range(len(pulse_amp)):
                if not np.isnan(pulse_amp[i]):
                    cntct_dicts[cntcts_made] = {
                        # Assuming one-indexed contact ids
                        "Contact_ID": i+1,
                        "Active": True,
                        pulse_amp_key: pulse_amp[i],
                        not_pulse_amp_key: False}
                    cntcts_made += 1
            
            elec_dicts[0] = {
                # Assuming both electrodes are the same type
                "Name": oss_elec_name,
                "PathToCustomParameters": "",
                "Rotation[Degrees]": self.get_rot_z(index_side),
                "Direction": {"x[mm]": unit_directions[index_side,0],
                              "y[mm]": unit_directions[index_side,1],
                              "z[mm]": unit_directions[index_side,2]},
                "TipPosition": {
                    "x[mm]": tip_pos[index_side,0],
                    "y[mm]": tip_pos[index_side,1],
                    "z[mm]": tip_pos[index_side,2]},
                "Contacts": cntct_dicts.tolist()}
          
            
        ### MAKE THE DICTIONARY ###    
        
        partial_dict = {
	    "ModelSide": 0,  # one per simulation
        "BrainRegion": {
        "Center":{"x[mm]":self.get_imp_coord()[hemis_idx,0],"y[mm]":self.get_imp_coord()[hemis_idx,1],"z[mm]":self.get_imp_coord()[hemis_idx,2]}},
        "Electrodes": elec_dicts.tolist(),
        # You want to activate a surface for grounding if there is a monopolar electrode in the model
        "Surfaces": [{
                "Name":"BrainSurface",
                "Active": bool(mono_elec),
                "Current[A]": 0.0,
                "Voltage[V]": 0.0}],
        "MaterialDistribution": {
            "MRIPath": self.get_mri_name(),
            "MRIMapping": {
                # Make UNK map to one lower than the indices provided in the lead settings
                UNK: int(min([self.get_csf_idx(),self.get_gm_idx(),self.get_wm_idx()])-1),
                CSF:self.get_csf_idx(),
                WM:self.get_wm_idx(),
                GM:self.get_gm_idx(),
                # Make BLOOD map to one higher than the indices provided in the lead settings
                BLOOD:int(max([self.get_csf_idx(),self.get_gm_idx(),self.get_wm_idx()])+1)},
            "DiffusionTensorActive": len(self.get_dti_name())>0,
            "DTIPath": self.get_dti_name()},
            "StimulationSignal": {"CurrentControlled": self.get_cur_ctrl()[hemis_idx]},
        "CalcAxonActivation": self.get_calc_axon_act(),
        "ActivationThresholdVTA": self.get_act_thresh_vta(),
        "OutputPath": os.path.join(os.environ["STIMFOLDER"],HEMIS_OUTPUT_PATHS[hemis_idx]),
        "TemplateSpace": self.get_est_in_temp()}


        partial_settings = Settings(partial_dict)
        return partial_settings.complete_settings()
    
    
    def save_to_oss_json(self, json_path, hemis_idx = 0):
        
        settings = self.make_settings(hemis_idx)
        with open(json_path, 'w') as f:
            json.dump(total_dict, f)
    
    def get_num_elecs(self):
        return self.NUM_ELECS
    
    def get_pat_fold(self):
        return self._get_str("Patient_folder")
    
    def get_est_in_temp(self):
        return bool(self._get_num("Estimate_In_Template"))
    
    # MRIPath
    def get_mri_name(self):
        return self._get_str("MRI_data_name")
    
    def get_dti_name(self):
        return self._get_str("DTI_data_name")
    
    def get_gm_idx(self):
        return int(self._get_num("GM_index"))
    
    def get_wm_idx(self):
        return int(self._get_num("WM_index"))
    
    def get_csf_idx(self):
        return int(self._get_num("CSF_index"))
    
    def get_def_mat(self):
        return self._get_str("default_material")
    
    def get_elec_type(self):
        return self._get_str("Electrode_type")
    
    def get_cntct_loc(self):
        e1 = np.asarray(self._file[self._settings["contactLocation"][0,0]][:,:])
        e2 = np.asarray(self._file[self._settings["contactLocation"][1,0]][:,:])
        return np.stack((e1,e2))
    
    # Used to re-compute rot_z 
    def get_y_mark_nat(self):
        return self._get_arr("yMarkerNative")
    
    def get_y_mark_mni(self):
        return self._get_arr("yMarkerMNI")
    
    def get_head_nat(self):
        return self._get_arr("headNative")
    
    def get_head_mni(self):
        return self._get_arr("headMNI")
    
    def get_imp_coord(self):
        return self._get_arr("Implantation_coordinate")
    
    def get_sec_coord(self):
        return self._get_arr("Second_coordinate")
    
    # Always recalculated from the other settings
    # IMPORTANT: it is actually not native but scrf!
    def get_rot_z(self,index_side):
        head_nat = self.get_head_nat()[:, index_side]
        y = self.get_y_mark_nat()[:, index_side] - head_nat
        y_postop = y / np.linalg.norm(y)
        phi = np.arctan2(-y_postop[0], y_postop[1])
        return phi * 180.0 / np.pi
    
    def get_stim_set_mode(self):
        return self._get_num("stimSetMode")
    
    def get_cur_ctrl(self):
        return self._get_arr("current_control").T[0]
    
    def get_phi_vec(self):
        return self._get_arr("Phi_vector")
    
    def get_case_grnd(self):
        return self._get_arr("Case_grounding")
    
    def get_act_thresh_vta(self):
        return self._get_num("Activation_threshold_VTA")
    
    def get_calc_axon_act(self):
        return self._get_num("calcAxonActivation")
    
    def get_connectome(self):
        return self._get_str("connectome")
    
    def get_axon_len(self):
        return self._get_arr("axonLength")
    
    def get_fib_diam(self):
            return self._get_arr("fiberDiameter")
    
    def get_conectome_path(self):
        return self._get_str("connectomePath")
    
    def get_connectome_tract_names(self):
        # Array of h5py references to the ascii arrays of the tract names
        ref_arr = self._settings['connectomeTractNames'][0]
        
        # Initialize the tract names array
        tract_names = np.empty(len(ref_arr),dtype=object)
   
        # Convert the references to name strings and return the strings
        for i in range(len(ref_arr)):
            entry = self._file[ref_arr[i]][:]
            ascii_codes = np.ndarray(entry.shape[0])
            for j in range(entry.shape[0]):
                ascii_codes[j] = entry[j][0]
            tract_names[i] = ''.join(np.vectorize(chr)(ascii_codes.astype(int)))
        return tract_names
    
    def get_inter_mode(self):
        return self._get_num("interactiveMode")
    
        
    # Private fxns
    
    def _is_h5(self):
        return self._h5
    
    # Worth considering making different subclasses for different file types?
    def _get_num(self,field_name):
        if self._is_h5():
            return self._settings[field_name][0][0]
        else:
            return self._settings[field_name][0][0][0][0]
    
    def _get_str(self,field_name):
        entry = ''
        if self._is_h5():
            ascii_codes = self._settings[field_name][:,0]
            entry = ''.join(np.vectorize(chr)(ascii_codes))
        else:
            entry = self._settings[field_name][0][0][0]
        if entry == 'no dti':
            entry = ''
        return entry
        
    def _get_arr(self,field_name):
        if self._is_h5():
            return self._settings[field_name][:,:].T
        else:
            if len(self._settings[field_name][0][0]) == 1:
                return self._settings[field_name][0][0][0]
            else:
                return self._settings[field_name][0][0].astype(float)