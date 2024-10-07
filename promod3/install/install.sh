PROMOD_VERSION="3.4.1"
SRC_FOLDER="../src"
OST_ROOT="/home/jerry/mara/protein_structure_processing/promod3/openstructure/stage"
mkdir $SRC_FOLDER

# ENVIRONMENT
##############################################################################
export PROMOD_VERSION="${PROMOD_VERSION}"
export PROMOD_ROOT="/home/jerry/mara/protein_structure_processing/promod3"


cd ${SRC_FOLDER} && \
# copy promod release
wget -O promod-${PROMOD_VERSION}.tar.gz -nc https://git.scicore.unibas.ch/schwede/ProMod3/-/archive/${PROMOD_VERSION}/promod-${PROMOD_VERSION}.tar.gz && \
mkdir promod-${PROMOD_VERSION} && \
tar xf promod-${PROMOD_VERSION}.tar.gz -C ${SRC_FOLDER}/promod-${PROMOD_VERSION} --strip-components=1 && \
mkdir -p ${SRC_FOLDER}/promod-${PROMOD_VERSION}/build && \
cd ${SRC_FOLDER}/promod-${PROMOD_VERSION}/build && \
# Build and install ProMod3
cmake .. -DOST_ROOT=${OST_ROOT} \
            -DOPTIMIZE=1 \
            -DENABLE_SSE=1 \
            -DDISABLE_DOCUMENTATION=1 && \

make -j4 && make check && make install && \
# cleanup
cd ${SRC_FOLDER} && rm ${SRC_FOLDER}/promod-${PROMOD_VERSION}.tar.gz && \
rm -rf ${SRC_FOLDER}/promod-${PROMOD_VERSION}
