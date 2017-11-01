using GenomicIndexes
using Base.Test

@testset "Tabix" begin
    Tabix = GenomicIndexes.Tabix
    tabixfile = joinpath(dirname(@__FILE__), "TAIR10_GFF3_genes.sorted.gff.bgz.tbi")

    # Read from a file.
    tbi = Tabix(tabixfile)
    @test tbi isa Tabix
    @test GenomicIndexes.seqnames(tbi) == ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM"]

    # Read from a stream.
    tbi = Tabix(open(tabixfile))
    @test tbi isa Tabix
    @test GenomicIndexes.seqnames(tbi) == ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM"]

    # Test overlapchunks.
    @test GenomicIndexes.overlapchunks(tbi, "Chr1", 240057:242608) == [
        (0x0000000000000000, 0x00000006dfe7e426),
        (0x0000000a9781af62, 0x0000000a9781b540)]
    @test GenomicIndexes.overlapchunks(tbi, "Chr4", 677054:679237) == [
        (0x0000002fa1fb670b, 0x0000003376fb1fbc)]
    @test GenomicIndexes.overlapchunks(tbi, "ChrM", 292973:293431) == [
        (0x0000004d38e83e1c, 0x0000004d38e8fa5c)]
    @test_throws ArgumentError GenomicIndexes.overlapchunks(tbi, "ChrZ", 1000:20000)
end

@testset "BAI" begin
    BAI = GenomicIndexes.BAI
    baifile = joinpath(dirname(@__FILE__), "ex1.bam.bai")

    bai = BAI(baifile)
    @test bai isa BAI

    bai = BAI(open(baifile))
    @test bai isa BAI
end

@testset "CSI" begin
    CSI = GenomicIndexes.CSI
    csifile = joinpath(dirname(@__FILE__), "TAIR10_GFF3_genes.sorted.gff.bgz.csi")

    # Load from a file.
    csi = CSI(csifile)
    @test csi isa CSI

    # Load from a stream.
    csi = CSI(open(csifile))
    @test csi isa CSI

    # Test overlapchunks.
    @test GenomicIndexes.overlapchunks(csi, 1, 240057:242608) == [
        (0x0000000000000000, 0x00000006dfe7e426),
        (0x0000000a9781af62, 0x0000000a9781b540)]
    @test GenomicIndexes.overlapchunks(csi, 4, 677054:679237) == [
        (0x0000002fa1fb670b, 0x0000003376fb1fbc)]
    @test GenomicIndexes.overlapchunks(csi, 7, 292973:293431) == [
        (0x0000004d38e83e1c, 0x0000004d38e8fa5c)]
    @test_throws ArgumentError GenomicIndexes.overlapchunks(csi, 10, 1000:20000)
end
