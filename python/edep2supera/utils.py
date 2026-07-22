import os
import edep2supera
from edep2supera import config
import ROOT
from ROOT import supera, TChain
supera.EDep
import numpy as np
from larcv import larcv

def dict2map(target):
    result = ROOT.std.map("std::string,std::string")()
    for key,val in target.items():
        result[str(key)] = str(val)
    return result

def get_edep2supera(config_key):

    e2s_driver = edep2supera.edep2supera.SuperaDriver()
    if os.path.isfile(config_key):
        e2s_driver.ConfigureFromFile(config_key)
    else:
        e2s_driver.ConfigureFromFile(config.get_config(config_key))
    
    return e2s_driver

def get_iomanager(outname):
    import tempfile
    cfg='''
IOManager: {
  Verbosity:   2
  Name:        "IOManager"
  IOMode:      1
  OutFileName: "%s" 
}
''' 
    #f=open('tmp.cfg','w')
    f=tempfile.NamedTemporaryFile(mode="w")
    f.write(cfg % outname)
    f.flush()
    o = larcv.IOManager(f.name)
    o.initialize()
    f.close()
    return o

def larcv_meta(supera_meta):
    larcv_meta = larcv.Voxel3DMeta()

    larcv_meta.set(supera_meta.min_x(),supera_meta.min_y(),supera_meta.min_z(),
                   supera_meta.max_x(),supera_meta.max_y(),supera_meta.max_z(),
                   supera_meta.num_voxel_x(),supera_meta.num_voxel_y(),supera_meta.num_voxel_z())
    
    return larcv_meta

def larcv_particle(p):
    
    US2NS = 1.e3

    larp=larcv.Particle()
    
    larp.id              (p.part.id)
    larp.shape           (int(p.part.shape))
    
    # particle's info setter
    larp.track_id         (p.part.trackid)
    if hasattr(p.part, "genid") and hasattr(larp, "gen_id"):
        larp.gen_id(p.part.genid)
        
    larp.pdg_code         (p.part.pdg)
    larp.momentum         (p.part.px,p.part.py,p.part.pz)
    larp.end_momentum     (p.part.end_px,p.part.end_py,p.part.end_pz)
    larp.distance_travel  (p.part.dist_travel)
    
    vtx_dict = dict(position = p.part.vtx, 
                    end_position = p.part.end_pt, 
                    first_step = p.part.first_step, 
                    last_step = p.part.last_step,
                    parent_position = p.part.parent_vtx,
                    ancestor_position = p.part.ancestor_vtx,
                   )
    for key,item in vtx_dict.items():
        getattr(larp,key)(item.pos.x, item.pos.y, item.pos.z, item.time * US2NS)
    
    #larp.distance_travel ( double dist ) { _dist_travel = dist; }
    larp.energy_init      (p.part.energy_init)
    larp.energy_deposit   (p.energy.sum())
    larp.creation_process (p.part.process)
    larp.num_voxels       (p.energy.size())
    
    # parent info setter
    larp.parent_track_id (p.part.parent_trackid)
    larp.parent_pdg_code (p.part.parent_pdg)
    larp.parent_creation_process(p.part.parent_process)
    larp.parent_id       (p.part.parent_id)
    for cid in p.part.children_id:
        larp.children_id(cid)

    # ancestor info setter
    larp.ancestor_track_id (p.part.ancestor_trackid)
    larp.ancestor_pdg_code (p.part.ancestor_pdg)
    larp.ancestor_creation_process(p.part.ancestor_process)
                                   
    if not p.part.group_id == supera.kINVALID_INSTANCEID: 
        larp.group_id(p.part.group_id)
        
    if not p.part.interaction_id == supera.kINVALID_INSTANCEID:
        larp.interaction_id(p.part.interaction_id)
    
    return larp

def run_supera(out_file='larcv.root',
    in_files=[],
    config_key='',
    num_events=-1,num_skip=0):
    
    reader = TChain("EDepSimEvents")
    in_files = list(in_files)
    print('Input edep-sim files:')
    for f in in_files:
        print(f)
        reader.AddFile(f)

    input_total = reader.GetEntries()
    if input_total < 1:
        print('No event found from',len(in_files),'files')
        for f in in_files:
            print(f)
        return 
    else:
        print('... processing',input_total,'events.')
        print('Output:',out_file)

    writer = get_iomanager(out_file)
    label_maker = get_edep2supera(config_key)
    
    id_vv=ROOT.std.vector("std::vector<unsigned long>")()
    value_vv=ROOT.std.vector("std::vector<float>")()

    id_v=ROOT.std.vector("unsigned long")()
    value_v=ROOT.std.vector("float")()
    
    if num_events < 0:
        num_events = reader.GetEntries()

    for e in range(reader.GetEntries()):

        if num_skip and e < num_skip:
            continue

        if num_events <= 0:
            break
        num_events -= 1
        
        reader.GetEntry(e)
        input_data = label_maker.ReadEvent(reader.Event)
        
        label_maker.GenerateImageMeta(input_data)
        label_maker.GenerateLabel(input_data)
        
        result = label_maker.Label()
        meta   = larcv_meta(label_maker.Meta())
        
        tensor_energy = writer.get_data("sparse3d","pcluster")
        result.FillTensorEnergy(id_v,value_v)
        larcv.as_event_sparse3d(tensor_energy,meta,id_v,value_v)

        tensor_semantic = writer.get_data("sparse3d","pcluster_semantics")
        result.FillTensorSemantic(id_v,value_v)
        larcv.as_event_sparse3d(tensor_semantic,meta,id_v,value_v)

        cluster_energy = writer.get_data("cluster3d","pcluster")
        result.FillClustersEnergy(id_vv,value_vv)
        larcv.as_event_cluster3d(cluster_energy,meta,id_vv,value_vv)

        cluster_dedx = writer.get_data("cluster3d","pcluster_dedx")
        result.FillClustersdEdX(id_vv,value_vv)
        larcv.as_event_cluster3d(cluster_dedx,meta,id_vv,value_vv)
        
        particle = writer.get_data("particle","pcluster")
        for p in result._particles:
            if not p.valid:
                continue
            larp = larcv_particle(p)
            particle.append(larp)

        writer.set_id(reader.Event.RunId,0,reader.Event.EventId)
        writer.save_entry()
        
    writer.finalize()
