class SimulationBox:
    def __init__(self,xlen=100, ylen=100, zlen=100, periodic="plane", custom={'x':'','y':'','z':''}):
        self.xlen = xlen
        self.ylen = ylen
        self.zlen = zlen
        self.periodic = periodic
        self.custom = custom

    def dim_string(self):
        return f'x:{self.xlen} y:{self.ylen} z:{self.zlen}'

    def volume(self):
        return (self.xlen*1e-6)*(self.ylen*1e-6)*(self.zlen*1e-6)