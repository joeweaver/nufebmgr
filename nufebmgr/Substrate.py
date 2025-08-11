from dataclasses import dataclass
from typing import Literal, Optional
@dataclass
class Substrate:
    name: str
    init_concentration: float
    bulk_concentration: float
    molecular_weight: Optional[float] = None