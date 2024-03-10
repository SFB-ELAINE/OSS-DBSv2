import fileinput
from abc import ABC, abstractmethod


class NeuronModel(ABC):
    """Interface to NEURON simulator."""

    @abstractmethod
    def paste_to_hoc(self):
        """Paste Python parameters into HOC file."""
        pass

    @abstractmethod
    def paste_paraview_vis(self):
        """Convert to paraview file."""


class MRG2002(NeuronModel):
    def paste_to_hoc(
        self,
        axonnodes,
        paranodes1,
        paranodes2,
        axoninter,
        axontotal,
        v_init,
        fiberD,
        paralength1,
        paralength2,
        nodelength,
        nodeD,
        axonD,
        paraD1,
        paraD2,
        deltax,
        nl,
        steps_per_ms,
    ):
        axonnodes_line = "axonnodes="
        axonnodes_input = f"axonnodes={axonnodes}\n"

        paranodes1_line = "paranodes1="
        paranodes1_input = f"paranodes1={paranodes1}\n"

        paranodes2_line = "paranodes2="
        paranodes2_input = f"paranodes2={paranodes2}\n"

        axoninter_line = "axoninter="
        axoninter_input = f"axoninter={axoninter}\n"

        axontotal_line = "axontotal="
        axontotal_input = f"axontotal={axontotal}\n"

        nv_init_line = "v_init="
        nv_init_input = f"v_init={v_init}\n"  # normally, -80mv

        fiberD_line = "fiberD="
        fiberD_input = f"fiberD={fiberD}\n"  # fiber diameter

        paralength1_line = "paralength1="
        paralength1_input = f"paralength1={paralength1}\n"

        paralength2_line = "paralength2="
        paralength2_input = f"paralength2={paralength2}\n"

        nodelength_line = "nodelength="
        nodelength_input = f"nodelength={nodelength}\n"

        nodeD_line = "nodeD="
        nodeD_input = f"nodeD={nodeD}\n"

        axonD_line = "axonD="
        axonD_input = f"axonD={axonD}\n"

        paraD1_line = "paraD1="
        paraD1_input = f"paraD1={paraD1}\n"

        paraD2_line = "paraD2="
        paraD2_input = f"paraD2={paraD2}\n"

        deltax_line = "deltax="
        deltax_input = f"deltax={deltax}\n"

        nl_line = "nl="
        nl_input = f"nl={nl}\n"

        steps_per_ms_line = "steps_per_ms="
        steps_per_ms_input = f"steps_per_ms={steps_per_ms}\n"

        x = fileinput.input(files="axon4pyfull.hoc", inplace=1)
        for line in x:
            if line.startswith(axonnodes_line):
                line = axonnodes_input
            if line.startswith(paranodes1_line):
                line = paranodes1_input
            if line.startswith(paranodes2_line):
                line = paranodes2_input
            if line.startswith(axoninter_line):
                line = axoninter_input
            if line.startswith(axontotal_line):
                line = axontotal_input
            if line.startswith(nv_init_line):
                line = nv_init_input
            if line.startswith(fiberD_line):
                line = fiberD_input
            if line.startswith(paralength1_line):
                line = paralength1_input
            if line.startswith(paralength2_line):
                line = paralength2_input
            if line.startswith(nodelength_line):
                line = nodelength_input

            if line.startswith(nodeD_line):
                line = nodeD_input
            if line.startswith(axonD_line):
                line = axonD_input
            if line.startswith(paraD1_line):
                line = paraD1_input
            if line.startswith(paraD2_line):
                line = paraD2_input
            if line.startswith(deltax_line):
                line = deltax_input
            if line.startswith(nl_line):
                line = nl_input
            if line.startswith(steps_per_ms_line):
                line = steps_per_ms_input
            print(line, end="")
        x.close()

        return True

    def paste_paraview_vis(self, Points_on_model, N_comp_in_between):
        NPoints_line = "Points_on_model"
        NPoints_input = f"Points_on_model={Points_on_model}\n"  # NEURON uses ms
        N_comp_in_between_line = "N_comp_in_between"
        N_comp_in_between_input = (
            f"N_comp_in_between={N_comp_in_between}\n"  # NEURON uses ms
        )

        fl = fileinput.input(
            files="Visualization_files/Paraview_vis_axon.py", inplace=1
        )
        for line in fl:
            if line.startswith(NPoints_line):
                line = NPoints_input
            if line.startswith(N_comp_in_between_line):
                line = N_comp_in_between_input
            print(line, end="")
        fl.close()

        return True


class McNeal1976(NeuronModel):
    def paste_to_hoc(self, axonnodes, axoninter, axontotal, v_init, steps_per_ms):
        NNODES_line = "NNODES ="
        NNODES_input = f"NNODES = {axonnodes}\n"

        axonnodes_line = "axonnodes="
        axonnodes_input = f"axonnodes={axonnodes}\n"

        nv_init_line = "v_init="
        nv_init_input = f"v_init={v_init}\n"

        steps_per_ms_line = "steps_per_ms="
        steps_per_ms_input = f"steps_per_ms={steps_per_ms}\n"

        x = fileinput.input(files="init_B5_extracellular.hoc", inplace=1)
        for line in x:
            if line.startswith(axonnodes_line):
                line = axonnodes_input
            if line.startswith(nv_init_line):
                line = nv_init_input
            if line.startswith(steps_per_ms_line):
                line = steps_per_ms_input
            print(line, end="")
        x.close()

        x = fileinput.input(files="axon5.hoc", inplace=1)
        for line in x:
            if line.startswith(NNODES_line):
                line = NNODES_input
            print(line, end="")
        x.close()

        return True

    def paste_paraview_vis(self, Points_on_model, N_comp_in_between):
        NPoints_line = "Points_on_model"
        NPoints_input = f"Points_on_model={Points_on_model}\n"  # NEURON uses ms
        N_comp_in_between_line = "N_comp_in_between"
        N_comp_in_between_input = (
            f"N_comp_in_between={N_comp_in_between}\n"  # NEURON uses ms
        )

        fl = fileinput.input(
            files="Visualization_files/Paraview_vis_axon.py", inplace=1
        )
        for line in fl:
            if line.startswith(NPoints_line):
                line = NPoints_input
            if line.startswith(N_comp_in_between_line):
                line = N_comp_in_between_input
            print(line, end="")
        fl.close()

        return True
