using GenomicIndexes
using Base.Test

@testset "Tabix" begin
    tabixfile = joinpath(dirname(@__FILE__), "TAIR10_GFF3_genes.sorted.gff.bgz.tbi")
    tabix = Tabix(tabixfile)
    @test tabix isa Tabix
end
