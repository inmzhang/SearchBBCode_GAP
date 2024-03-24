# Install HPC-GAP
install-gap:
    ./scripts/install_hpc_gap.sh

# Search BB codes with GAP
# search num_threads='8':
#     time gap -P {{num_threads}} --quitonbreak -b -q search.g -c 'QUIT;'