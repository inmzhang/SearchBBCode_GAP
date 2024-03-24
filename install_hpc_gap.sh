# Install the HPC-GAP system(dev-build)
set -e

git clone https://github.com/gap-system/gap hpcgap

cd hpcgap

./autogen.sh
./configure --enable-hpcgap
make

make bootstrap-pkg-full

# tomlib is broken in the dev-build
# rm -rf pkg/tomlib

sudo ln -s $PWD/gap /usr/local/bin
