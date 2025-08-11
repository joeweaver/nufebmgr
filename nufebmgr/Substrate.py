from dataclasses import dataclass
from typing import Literal, Optional
@dataclass
class Substrate:
    name: str
    init_concentration: float
    bulk_concentration: float
    diffusion_coefficient: float
    biofilm_diffusion_ratio: float
    molecular_weight: Optional[float] = None