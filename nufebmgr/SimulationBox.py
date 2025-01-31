class SimulationBox:
    def __init__(self,xlen=100, ylen=100, zlen=100, periodic="plane", custom={'x':'','y':'','z':''}):
        self.xlen = xlen
        self.ylen = ylen
        self.zlen = zlen
        self.periodic = periodic
        self.custom = custom
