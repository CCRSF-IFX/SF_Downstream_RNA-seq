my_envs = ['R.yaml','multiqc.yaml', 'myenv.yaml']

rule make_all_envs:
    input:
        expand("created-{name}", name=my_envs)


for env_file in my_envs:
    rule:
        output:
            temp("created-%s" % env_file)
        conda:
            "envs/%s" % env_file
        shell:
            "touch {output}"
