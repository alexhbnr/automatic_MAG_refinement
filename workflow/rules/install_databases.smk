#### Install databases #########################################################

rule install_dbs:
    input:
        mmseqs2_gtdb = "{resourcesdir}/mmseqs2/gtdb/mmseqs2_gtdb_mapping",
        checkm_db = "{resourcesdir}/checkM/setRoot.done",
        gunc_db = "{resourcesdir}/GUNC/db"
    output:
        touch("{resourcesdir}/databases_installation.done")

rule install_mmseqs2_gtdb_db:
    output:
        "{resourcesdir}/mmseqs2/gtdb/mmseqs2_gtdb_mapping"
    message: "Install GTDB database for MMSeqs2 aminoacid profiling"
    conda: "envs/MMSeqs2_GTDB.yaml"
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
    conda: "envs/checkM.yaml"
    params:
        dbdir = "{resourcesdir}/checkM"
    shell:
        """
        echo {params.dbdir} | checkm data setRoot {params.dbdir}
        """

rule install_gunc_database:
    output:
        directory("{resourcesdir}/GUNC/db")
    message: "Install GUNC database"
    conda: "envs/GUNC.yaml"
    params:
        dir = "{resourcesdir}/GUNC/db"
    shell:
        """
        mkdir {params.dir}
        gunc download_db {params.dir}
        """

################################################################################

