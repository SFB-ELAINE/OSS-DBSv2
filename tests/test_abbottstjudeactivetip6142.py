from ossdbs.electrodes.abbott_stjude_active_tip_6142_6145 import AbbottStjudeActiveTip6142_6145
import pytest
from tests.electrode_testbase import ElectrodeTestBase


@pytest.mark.skip
class TestAbbottStJudeActiveTip6142_6145(ElectrodeTestBase):
    FILE_PREFIX="AbbottStjudeActiveTip6142_6145"
    TESTDATA = [
        ((0, 0, 0), (0, 0, 1), 0.0,FILE_PREFIX+ '_1.json'),
        ((0, 0, 0), (0, 0, 0), 0.0,FILE_PREFIX+ '_1.json'),
        ((0, 0, 0), (0, 0, 1), 3.0,FILE_PREFIX+ '_1.json'),
        ((1, -2, 3), (0, 0, 1), 0.0,FILE_PREFIX+ '_2.json'),
        ((1, -2, 3), (0, 0, 0), 0.0,FILE_PREFIX+ '_2.json'),
        ((1, -2, 3), (0, 0, 0), 3.0,FILE_PREFIX+ '_2.json'),
        ((1, -2, 3), (2.0, 0, 1.0), 0.0,FILE_PREFIX+'_3.json'), 
        # Behaviour depends on python handling the floating point operations
        ((1, -2, 3), (2.0/3.0, 0, 1.0/3.0), 0.0,FILE_PREFIX+'_3.json'),
        ((1, -2, 3), (2.0/3.0, 0, 1.0/3.0), 3.0,FILE_PREFIX+'_3.json'),
    ]

    @pytest.mark.parametrize('translation, direction, rotation, path',
                             TESTDATA)
    def test_generate_geometry(self,
                               translation,
                               direction,
                               rotation,
                               path) -> None:
        electrode = AbbottStjudeActiveTip6142_6145(translation,
                                                   direction,
                                                   rotation)
        desired = self.load_json(path=path)
        geometry = electrode.generate_geometry()
        assert desired == self.electrode_to_dict(geometry)

    def test_generate_geometry_default(self):
        electrode = AbbottStjudeActiveTip6142_6145()
        desired = self.load_json(path='case1.json')
        geometry = electrode.generate_geometry()
        assert desired == self.electrode_to_dict(geometry)
