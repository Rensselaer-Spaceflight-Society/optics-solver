
# STAGE CLASS
# okay so the stage essentially acts as a sort of container for all the optics
# run a command to render and evaluate all the outputs of the "diagnostic" elements
# possibly add a child class for transients
# - Possible "OptimizationStage" child class that lets you optimize a setup for a certain goal
# ^ This is where something like element overlap might matter
# Another Stage note -> if you use a transient stage, you'd could use a minkowski metric with rotations?
import copy
from models.diagnostic import Diagnostic

class Stage:
    AXIS = 2  # which coordinate of position.coords is "along the arm"

    def __init__(self, arms=0, dispersion=None):
        self.arms = {i: {} for i in range(arms)}

    def set_arm(self, id, params):
        self.arms[id] = {
            "max_length": params["length"],
            "source": params["source"],   # template source, never mutated
            "objects": [],
            "misc_params": {k: v for k, v in params.items() if k not in ("length", "source")},
        }

    def add_obj(self, arm_id, obj):
        self.arms[arm_id]["objects"].append(obj)

    def run(self):
        results = {}
        for arm_id, arm in self.arms.items():
            src = copy.deepcopy(arm["source"])# transfer() changes the source
            objs = sorted(arm["objects"], key=lambda o: o.position.coords[self.AXIS])
            for o in objs:
                if isinstance(o, Diagnostic):
                    o.data.clear()# reset readings/stuff

            z = 0.0
            for o in objs:
                oz = o.position.coords[self.AXIS]
                if oz < z:
                    raise ValueError(f"Arm {arm_id}: element at z={oz} overlaps the previous one (ends at z={z})")
                if oz > arm["max_length"]:
                    raise ValueError(f"Arm {arm_id}: element at z={oz} is past max_length")
                # gap since the last element
                if oz > z:
                    src.apply(src.translation(oz - z, src.n))
                src = o.transfer(src)
                z = oz + getattr(o, "thickness", 0.0)# lens matrix already includes its thickness so dw abt it
            results[arm_id] = {o: o.data[-1] for o in objs if isinstance(o, Diagnostic)}
        return results

