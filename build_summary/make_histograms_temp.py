import os
import glob
import pr500_qc as pr
dep_path = '/Users/elemire/Workspace/deprecated/devwork'

dets = list(set(os.path.basename(f)[0:14] for f in glob.glob(os.path.join('/Users/elemire/Workspace/deprecated/devwork', 'P*'))))
for det in dets:
    print det
    pr.median_viability_histogram(os.path.join(dep_path, det + '_NORM.gct'),
                                  os.path.join(dep_path, 'QC', det, '{}_median_viability_dist.png'.format(det)), True)



import sys
sys.path.append('/Users/elemire/Workspace/prism_pipeline')
import os
import glob
import pr500_qc as pr

dep_path = '/Users/elemire/Workspace/deprecated/devwork'
dets = list(set(os.path.basename(f)[0:14] for f in glob.glob(os.path.join('/Users/elemire/Workspace/deprecated/devwork', 'P*'))))

for det in dets:
    pr.distributions(dep_path, [det], os.path.join(dep_path, 'QC', det))






pr.distributions('GCT', ['PCAL102_CS5_X1', 'PCAL102_CS5_X2', 'PCAL102_CS5_X3', 'PCAL102_CS5_X4', 'PCAL102_CS5_X5', 'PCAL102_CS5_X6'], 'GCT/QC')