Bootstrap: docker
From: conda/miniconda3

%files
	environment.yml /opt/environment.yml

%post
	apt-get update && apt-get install build-essential -y
	
	conda env create -f /opt/environment.yml && conda clean -afy

	export PATH="/usr/local/envs/pagoo/bin:$PATH"

	R --slave -e 'devtools::install_github("iferres/pagoo")'
	R --slave -e 'install.packages("rhierbaps")'

%environment
	export PATH="/usr/local/envs/pagoo/bin:$PATH"

%labels
	Authors: Ignacio Ferres (iferres@pasteur.edu.uy) & Gregorio Iraola (giraola@pasteur.edu.uy)
	Maintainer: Ignacio Ferres (iferres@pasteur.edu.uy)
