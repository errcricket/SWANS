#config files
configfile: 'configs/local_configs.yaml' 

#paths

cellmarker_dir = 'reference_data/CellMarker/'
panglaoDB_dir = 'reference_data/panglaoDB/'

#cellmarkerDB marker files
cellmarker_config = {'human':cellmarker_dir + 'cell_types/human/', 'mouse':cellmarker_dir + 'cell_types/mouse/'}
cellmarker_config[ config['ORGANISM'] ]
CELLMARKER_CONFIG = cellmarker_config[ config['ORGANISM'] ]

#cellmarker gmt file
cellmarker_gmt_config = {'human':cellmarker_dir + 'CellMarker_human.gmt', 'mouse':cellmarker_dir + 'CellMarker_mouse.gmt'}
cellmarker_gmt_config[ config['ORGANISM'] ]
CELLMARKER_GMT_CONFIG = cellmarker_gmt_config[ config['ORGANISM'] ]

#panglaoDB marker files
panglaoDB_config = {'human':panglaoDB_dir + 'cell_types/human/', 'mouse':panglaoDB_dir + 'cell_types/mouse/'}
panglaoDB_config[ config['ORGANISM'] ]
PANGLAODB_CONFIG = panglaoDB_config[ config['ORGANISM'] ]

#panglaoDB gmt file
panglaoDB_gmt_config = {'human':panglaoDB_dir + 'panglaoDB_human.gmt', 'mouse':panglaoDB_dir + 'panglaoDB_mouse.gmt'}
panglaoDB_gmt_config[ config['ORGANISM'] ]
PANGLAODB_GMT_CONFIG = panglaoDB_gmt_config[ config['ORGANISM'] ]
