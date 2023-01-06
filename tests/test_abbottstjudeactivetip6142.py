

from ossdbs.electrodes.abbott_stjude_active_tip_6142_6145 import AbbottStjudeActiveTip6142_6145
import pytest


@pytest.mark.skip
class TestAbbottStJudeeActiveTip6142_6145:

    TESTDATA = [
        ((0, 0, 0), (0, 0, 1), 0.0, 'case1.json'),
        ((0, 0, 0), (0, 0, 0), 0.0, 'case1.json'),
        ((0, 0, 0), (0, 0, 1), 3.0, 'case1.json'),
        ((1, -2, 3), (0, 0, 1), 0.0, 'case2.json'),
        ((1, -2, 3), (0, 0, 1), 0.0, 'case3.json'),
    ]

    @pytest.mark.parametrize('translation, direction, rotation, path',
                             TESTDATA)
    def test_generate_geometry(self,
                               translation,
                               direction,
                               rotation,
                               path):
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
