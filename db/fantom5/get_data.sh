wget http://enhancer.binf.ku.dk/presets/facet_expressed_enhancers.tgz
mkdir organ_enhancers
tar -xzf facet_expressed_enhancers.tgz -C organ_enhancers/ --wildcards UBERON*
# celltype_enhancers
# tar -xzf facet_expressed_enhancers.tgz -C celltype_enhancers/ --wildcards CL:*

