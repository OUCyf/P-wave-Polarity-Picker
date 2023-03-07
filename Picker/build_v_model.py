
from obspy.taup.taup_create import build_taup_model

v_model_path  = './v_model.nd'
build_taup_model(v_model_path, output_folder='./', verbose=True)
