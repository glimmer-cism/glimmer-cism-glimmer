PYTHON=python
GLIDE=simple_glide
TEST_OUTPUT=ishom.a.out.nc ishom.b.out.nc ishom.c.out.nc ishom.d.out.nc

all: ${TEST_OUTPUT}

ishom.a.nc: ishom.a.config ismip_hom_a.py glimcdf.py
	${PYTHON} ismip_hom_a.py

ishom.a.out.nc: ishom.a.config ishom.a.nc 
	echo ishom.a.config | ${GLIDE}

ishom.b.nc: ishom.b.config ismip_hom_b.py glimcdf.py
	${PYTHON} ismip_hom_b.py

ishom.b.out.nc: ishom.b.config ishom.b.nc
	echo ishom.b.config | ${GLIDE}

ishom.c.nc: ishom.c.config ismip_hom_c.py glimcdf.py
	${PYTHON} ismip_hom_c.py

ishom.c.out.nc: ishom.c.config ishom.c.nc
	echo ishom.c.config | ${GLIDE}

ishom.d.nc: ishom.d.config ismip_hom_d.py glimcdf.py
	${PYTHON} ismip_hom_d.py

ishom.d.out.nc: ishom.d.config ishom.d.nc
	echo ishom.d.config | ${GLIDE}

ishom.f.nc: ishom.f.config ismip_hom_f.py glimcdf.py
	${PYTHON} ismip_hom_f.py

ishom.f.out.nc: ishom.f.config ishom.f.nc
	echo ishom.f.config | ${GLIDE}


