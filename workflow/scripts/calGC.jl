using FASTX
using BioSequences

function calGC(seq::LongSequence{DNAAlphabet{4}})::Float64
    BioSequences.count_naive(x->x==DNA_G || x==DNA_C,  seq)/length(seq) *100
end

function main(fa_fname::String,bedFile::String,outFile::String)
    fasta = FASTA.Reader(open(fa_fname),index=string(fa_fname,".fai"))
    open(outFile,"w") do l
        println(l,"chr\tstart\tstop\tgc")
        for line in eachline(bedFile)
            if startswith(line,"chrom")
                continue
            else
                chr,start,stop,_... = split(line,"\t")
                start = parse(Int64,start)
                stop = parse(Int64,stop)
                chr = convert(String,chr)
                seq = LongSequence{DNAAlphabet{4}}(extract(fasta,chr,start:stop))
                gc = calGC(seq)
                println(l,gc)
            end
        end
    end
    close(fasta)
end

main(snakemake.config["ref"],snakemake.input["ecc"],snakemake.output[1])
