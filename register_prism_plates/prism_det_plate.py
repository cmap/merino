

class PrismDetPlate(object):
    def __init__(self, pert_plate=None, davepool_id=None):
        self.pert_plate = pert_plate
        self.davepool_id = davepool_id
        self.pool_ids = set()
        self.rep_num = 0