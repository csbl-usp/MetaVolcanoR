context("ggplot2 objects as outputs of MetaVolcanoR")

data(diffexplist)

test_that("MetaVolcano outputs are MetaVolcano valid objects", {
    
    mv <- votecount_mv(diffexplist)
    expect_is(mv, "MetaVolcano")
    expect_is(mv@MetaVolcano, "gg")

    mv <- combining_mv(diffexplist)
    expect_is(mv, "MetaVolcano")
    expect_is(mv@MetaVolcano, "gg")

    mv <- rem_mv(diffexplist)
    expect_is(mv, "MetaVolcano")
    expect_is(mv@MetaVolcano, "gg")

})
