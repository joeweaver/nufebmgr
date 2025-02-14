import pytest
from nufebmgr.SimulationBox import SimulationBox

def test_volume_calc():
    s = SimulationBox() # 100 x 100 x 100 default
    assert s.volume() == pytest.approx(1e-12,abs=1e-14)

    s = SimulationBox(xlen=20,ylen=200,zlen=45)
    assert s.volume() == pytest.approx(1.8e-13, abs=1e-14)





