# this can definitely be cleaned up with wildcards!
# note that this file is written to work with the OSX 64-bit compiler!

APP1 = ephys_pipeline_cluster_standalone
APP2 = ephys_pipeline_lfp_standalone
APP3 = ephys_pipeline_mua_standalone
APP4 = ephys_pipeline_smscore_standalone
APP5 = ephys_pipeline_songextract_aggregate
APP6 = ephys_pipeline_sua_postproc_standalone
APP7 = ephys_pipeline_sua_standalone
APP8 = ephys_pipeline_sua_track_standalone

DIR = compiled
FIXPATH = /Contents/MacOS/

all : binary

binary : $(DIR)/$(APP1)/$(APP1).app/$(FIXPATH)/$(APP1) \
	$(DIR)/$(APP2)/$(APP2).app/$(FIXPATH)/$(APP2) \
	$(DIR)/$(APP3)/$(APP3).app/$(FIXPATH)/$(APP3) \
	$(DIR)/$(APP4)/$(APP4).app/$(FIXPATH)/$(APP4) \
	$(DIR)/$(APP5)/$(APP5).app/$(FIXPATH)/$(APP5) \
	$(DIR)/$(APP6)/$(APP6).app/$(FIXPATH)/$(APP6) \
	$(DIR)/$(APP7)/$(APP7).app/$(FIXPATH)/$(APP7) \
	$(DIR)/$(APP8)/$(APP8).app/$(FIXPATH)/$(APP8) \

$(DIR)/$(APP1)/$(APP1).app/$(FIXPATH)/$(APP1) : $(APP1).m
	mkdir -p $(DIR)/$(APP1)
	matlab -nodesktop -nosplash -r 'try,mcc -vm $(APP1).m -d $(DIR)/$(APP1);catch,end,exit'

$(DIR)/$(APP2)/$(APP2).app/$(FIXPATH)/$(APP2) : $(APP2).m ../ephys_visual_lfp_amp.m ../ephys_visual_lfp_tf.m
	mkdir -p $(DIR)/$(APP2)
	matlab -nodesktop -nosplash -r 'try,mcc -vm $(APP2).m -d $(DIR)/$(APP2);catch,end,exit'

$(DIR)/$(APP3)/$(APP3).app/$(FIXPATH)/$(APP3) : $(APP3).m ../ephys_visual_mua.m
	mkdir -p $(DIR)/$(APP3)
	matlab -nodesktop -nosplash -r 'try,mcc -vm $(APP3).m -d $(DIR)/$(APP3);catch,end,exit'

$(DIR)/$(APP4)/$(APP4).app/$(FIXPATH)/$(APP4) : $(APP4).m
	mkdir -p $(DIR)/$(APP4)
	matlab -nodesktop -nosplash -r 'try,mcc -vm $(APP4).m -d $(DIR)/$(APP4);catch,end,exit'

$(DIR)/$(APP5)/$(APP5).app/$(FIXPATH)/$(APP5) : $(APP5).m
	mkdir -p $(DIR)/$(APP5)
	matlab -nodesktop -nosplash -r 'try,mcc -vm $(APP5).m -d $(DIR)/$(APP5);catch,end,exit'

$(DIR)/$(APP6)/$(APP6).app/$(FIXPATH)/$(APP6) : $(APP6).m
	mkdir -p $(DIR)/$(APP6)
	matlab -nodesktop -nosplash -r 'try,mcc -vm $(APP6).m -d $(DIR)/$(APP6);catch,end,exit'

$(DIR)/$(APP7)/$(APP7).app/$(FIXPATH)/$(APP7) : $(APP7).m ../ephys_visual_sua.m ../helpers/ephys/ephys_spike_cluster_auto.m ../helpers/ephys/ephys_spike_clustergui_tetrode.m \
	../helpers/ephys/ephys_spike_detect.m ../helpers/ephys/get_geometric_coefficients.m ../helpers/ephys/get_wavelet_coefficients.m
	mkdir -p $(DIR)/$(APP7)
	matlab -nodesktop -nosplash -r 'try,mcc -vm $(APP7).m -d $(DIR)/$(APP7);catch,end,exit'

$(DIR)/$(APP8)/$(APP8).app/$(FIXPATH)/$(APP8) : $(APP8).m
	mkdir -p $(DIR)/$(APP8)
	matlab -nodesktop -nosplash -r 'try,mcc -vm $(APP8).m -d $(DIR)/$(APP8);catch,end,exit'




