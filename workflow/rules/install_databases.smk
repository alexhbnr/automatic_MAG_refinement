#### Install databases #########################################################

rule install_dbs:
    input:
        mmseqs2_gtdb = "{resourcesdir}/mmseqs2/gtdb/mmseqs2_gtdb_r207_db_mapping",
        checkm_db = "{resourcesdir}/checkM/setRoot.done",
        gunc_db = "{resourcesdir}/GUNC/db",
        checkm2_db = "{resourcesdir}/checkm2/CheckM2_database/uniref100.KO.1.dmnd",
        bakta = "{resourcesdir}/bakta/downloaded"
    output:
        touch("{resourcesdir}/databases_installation.done")

rule install_mmseqs2_gtdb_db:
    output:
        "{resourcesdir}/mmseqs2/gtdb/mmseqs2_gtdb_r207_db_mapping"
    message: "Install GTDB database for MMSeqs2 aminoacid profiling"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        cores = 32,
        mem = 32
    params: 
        outdir = "{resourcesdir}/mmseqs2/gtdb",
        tmpdir = config['tmpdir']
    threads: 32
    shell:
        """
        mmseqs databases GTDB  \
            {params.outdir} \
            {params.tmpdir} \
            --threads {threads}
        """

rule checkM_prepareDatabase:
    output:
        "{resourcesdir}/checkM/.dmanifest"
    message: "Download the checkM databases and extract the tarballs"
    params:
        outdir = "{resourcesdir}/checkM",
        url = "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
    shell:
        """
        cd {params.outdir}
        wget {params.url}
        tar xvf checkm_data_2015_01_16.tar.gz
        """

rule checkM_setRoot:
    input:
        "{resourcesdir}/checkM/.dmanifest"
    output:
        touch("{resourcesdir}/checkM/setRoot.done")
    message: "Specify the database folder in checkM"
    container: "https://depot.galaxyproject.org/singularity/checkm-genome:1.2.2--pyhdfd78af_1"
    params:
        dbdir = "{resourcesdir}/checkM"
    shell:
        """
        echo {params.dbdir} | checkm data setRoot {params.dbdir}
        """

rule install_checkm2_database:
    output:
        "{resourcesdir}/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
    message: "Install CheckM2 database"
    container: "https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0"
    params:
        dir = "{resourcesdir}/checkm2"
    shell:
        "checkm2 database --download --path {params.dir}"

rule install_gunc_database:
    output:
        directory("{resourcesdir}/GUNC/db")
    message: "Install GUNC database"
    container: "https://depot.galaxyproject.org/singularity/gunc:1.0.5--pyhdfd78af_0"
    params:
        dir = "{resourcesdir}/GUNC/db"
    shell:
        """
        mkdir {params.dir}
        gunc download_db {params.dir}
        """

rule download_bakta_db:
    output:
        touch("{resourcesdir}/bakta/downloaded")
    message: "Download and prepare the reference DB for Bakta"
    container: "/mnt/archgen/tools/singularity/containers/depot.galaxyproject.org-singularity-bakta-1.7.0--pyhdfd78af_0.img"
    resources:
        mem = 4,
        cores = 1
    params:
        dir = "{resourcesdir}/bakta"
    shell:
        """
        bakta_db download --output {params.dir}
        """

rule gtdbtk_download_db:
    output:
        f"{config['resourcesdir']}/gtdbtk/gtdbtk_{config['gtdb_version']}/metadata/metadata.txt"
    message: "Download and set-up the GTDBTK database"
    container: "docker://quay.io/biocontainers/gtdbtk:2.3.2--pyhdfd78af_0"
    params:
        url = {'r207_v2': "https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz",
               'r214.1': "https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz"}[config['gtdb_version']],
        tarball = f"{config['resourcesdir']}/gtdbtk/gtdbtk_{config['gtdb_version']}_data.tar.gz",
        gtdbtk_dir = "{config['resourcesdir']}/gtdbtk/gtdbtk_{config['gtdb_version']}"
    shell:
        """
        wget -O {params.tarball} {params.url} && \
        tar -xvzf {params.tarball} -C '{params.gtdbtk_dir}' --strip 1 > /dev/null && \
        rm {params.tarball} &&  \
        conda env config vars set GTDBTK_DATA_PATH='{gtdbtk_dir}'
        """

################################################################################

