String jvarkit(name) {
	return "java  -Djava.io.tmpdir=TMP -jar ${params.jvarkit} ${name} ";
	}
workflow {
	rotavirus_bams = Channel.of(
		"${params.testDir}/S1.bam",
		"${params.testDir}/S2.bam",
		"${params.testDir}/S3.bam",
		"${params.testDir}/S4.bam",
		"${params.testDir}/S5.bam"
		)
	
	rotavirus_vcfs = Channel.of(
		"${params.testDir}/rotavirus_rf.vcf.gz"
		).mix(rotavirus_bams.map{it.replaceAll("\\.bam",".vcf.gz")})

	TEST_BAMLEFTALIGN("${params.testDir}/rotavirus_rf.fa",rotavirus_bams.combine(Channel.of("--filter none","--filter only ","--filter discard","--regions RF02:1-10")))
	TEST_DICT2BED(rotavirus_bams.combine(Channel.of("")))
	TEST_FINDALLCOVERAGEATPOS("${params.testDir}/rotavirus_rf.fa",rotavirus_bams.collect(),Channel.of("-p 'RF01:100'"))
	TEST_GTF2BED(Channel.fromPath("${params.testDir}/Homo_sapiens.GRCh37.87.gtf.gz").combine(Channel.of("-c 'return gene_id' ")))
	TEST_SAMJDK("${params.testDir}/rotavirus_rf.fa",rotavirus_bams.combine(Channel.of("-e 'return record.getStart()<100;' ")))
	TEST_SAM2TSV("${params.testDir}/rotavirus_rf.fa",rotavirus_bams)
	TEST_VCFHEAD(rotavirus_vcfs.combine(Channel.of("-n 0","-n 10")))
	TEST_VCFFILTERSO(Channel.fromPath("${params.testDir}/rotavirus_rf.ann.vcf.gz").combine(Channel.of("-A 'SO:0001818,SO:0001629'  ")))
	TEST_VCFPOLYX("${params.testDir}/rotavirus_rf.fa",rotavirus_vcfs)
	TEST_VCFTAIL(rotavirus_vcfs.combine(Channel.of("-n 0","-n 10")))
	TEST_VCF2TABLE(rotavirus_vcfs.combine(Channel.of("--hide 'HOM_REF' ")))
	}


process TEST_BAMLEFTALIGN {
	tag "${bam}"
	afterScript "rm -rf TMP"
	input:
		val(fasta)
		tuple val(bam),val(args)
	script:
	"""
	mkdir -p TMP
	${jvarkit("bamleftalign")} ${args} -R '${fasta}' "${bam}" > TMP/jeter.sam
	"""
	}

process TEST_SAMJDK {
	tag "${bam}"
	afterScript "rm -rf TMP"
	input:
		val(fasta)
		tuple val(bam),val(args)
	script:
	"""
	mkdir -p TMP
	${jvarkit("samjdk")} ${args} -R '${fasta}' "${bam}" > TMP/jeter.sam
	"""
	}

process TEST_DICT2BED {
	tag "${bam}"
	afterScript "rm -rf TMP"
	input:
		tuple val(bam),val(args)
	script:
	"""
	mkdir -p TMP
	${jvarkit("dict2bed")} ${args} "${bam}" > TMP/jeter.bed
	"""
	}

process TEST_FINDALLCOVERAGEATPOS {
	input:
		val(genome)
		val(L)
		val(args)
	script:
	"""
	mkdir -p TMP
cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF
	cat TMP/jeter.list |\\
		 ${jvarkit("findallcoverageatposition")} -R ${genome} ${args} > TMP/jeter.tsv
	"""
	}

process TEST_GTF2BED {
        tag "${gtf}"
        afterScript "rm -rf TMP"
        input:
                tuple path(gtf),val(args)
        script:
        """
        mkdir -p TMP
        ${jvarkit("gtf2bed")} ${args} "${gtf}" > TMP/jeter.bed
        """
        }


process TEST_SAM2TSV {
	tag "${bam}"
	afterScript "rm -rf TMP"
	input:
		val(fasta)
		val(bam)
	script:
	"""
	mkdir -p TMP
	${jvarkit("sam2tsv")} -R '${fasta}' "${bam}" > TMP/jeter.bed
	"""
	}


process TEST_VCFHEAD {
	tag "${vcf}"
	afterScript "rm -rf TMP"
	input:
		tuple path(vcf),val(args)
	script:
	"""
	mkdir -p TMP
	${jvarkit("vcfhead")} ${args} "${vcf}" > TMP/jeter.vcf
	"""
	}

process TEST_VCFTAIL {
	tag "${vcf}"
	afterScript "rm -rf TMP"
	input:
		tuple path(vcf),val(args)
	script:
	"""
	mkdir -p TMP
	${jvarkit("vcftail")} ${args} "${vcf}" > TMP/jeter.vcf
	"""
	}


process TEST_VCFFILTERSO {
	tag "${vcf}"
	afterScript "rm -rf TMP"
	input:
		tuple path(vcf),val(args)
	script:
	"""
	mkdir -p TMP
	${jvarkit("vcffilterso")} ${args} "${vcf}" > TMP/jeter.vcf
	"""
	}

process TEST_VCFPOLYX {
	tag "${vcf}"
	afterScript "rm -rf TMP"
	input:
		val(fasta)
		val(vcf)
	script:
	"""
	mkdir -p TMP
	${jvarkit("vcfpolyx")} -n 2 -R '${fasta}' "${vcf}" > TMP/jeter.vcf
	"""
	}


process TEST_VCF2TABLE {
	tag "${vcf}"
	afterScript "rm -rf TMP"
	input:
		tuple val(vcf),val(args)
	script:
	"""
	mkdir -p TMP
	${jvarkit("vcf2table")} ${args} "${vcf}" > TMP/jeter.txt
	"""
	}
