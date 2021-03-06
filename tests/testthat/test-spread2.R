test_that("spread2 tests", {
  library(raster)
  library(data.table)
  library(fpCompare)
  library(quickPlot)
  on.exit(detach("package:raster"), add = TRUE)
  on.exit(detach("package:data.table"), add = TRUE)
  on.exit(detach("package:fpCompare"), add = TRUE)

  # inputs for x
  a <- raster(extent(0, 10, 0, 10), res = 1)
  b <- raster(a)
  sp <- 0.225
  spRas <- gaussMap(b)
  spRas[] <- spRas[] / maxValue(spRas) * sp / 2 + sp / 2 * 1.5
  b[] <- 1
  bb <- focal(b, matrix(1 / 9, nrow = 3, ncol = 3), fun = sum, pad = TRUE, padValue = 0)
  innerCells <- which(bb[] %==% 1)
  sams <- sample(innerCells, 2)

  set.seed(123)
  for (i in 1:20) {
    sams <- sample(innerCells, 2)
    out <- spread2(a, start = sams, spreadProb = 0.225, asRaster = FALSE)
    expect_true(length(unique(out$initialPixels)) == 2)
    expect_true(all(out$active == 0))
  }

  if (interactive()) print("testing maxSize")
  maxSizes <- 2:3
  for (i in 1:20) {
    seed <- sample(1e6, 1)
    set.seed(seed)
    sams <- sample(innerCells, 2)
    out <- spread2(a, start = sams, spreadProb = 0.225, maxSize = maxSizes,
                   asRaster = FALSE)
    expect_true(all(out[, .N, by = "initialPixels"]$N <= maxSizes[order(sams)]))
  }

  if (interactive()) print("testing exactSize")
  exactSizes <- c(5, 3.1)
  for (i in 1:20) {
    sams <- sample(innerCells, 2)
    out <- spread2(
      a,
      start = sams,
      spreadProb = 0.225,
      exactSize = exactSizes,
      asRaster = FALSE
    )
    attrib <- attr(out, "spreadState")$cluster$numRetries > 10
    if (any(attrib)) {
      frequ <- out[, .N, by = "initialPixels"]$N
      expect_true(all(frequ[attrib] <= floor(exactSizes[order(sams)][attrib])))
      expect_true(all(frequ[!attrib] == floor(exactSizes[order(sams)][!attrib])))
    } else {
      expect_true(all(out[, .N, by = "initialPixels"]$N == floor(exactSizes[order(sams)])))
    }
  }

  if (interactive()) print("testing exactSize")
  exactSizes <- c(5.01, 3.1, 4)
  for (i in 1:20) {
    sams <- sample(innerCells, length(exactSizes))
    out <- spread2(a, start = sams, spreadProb = 0.225, exactSize = exactSizes, asRaster = FALSE)
    attrib <- attr(out, "spreadState")$clusterDT$numRetries > 10
    if (any(attrib)) {
      frequ <- out[, .N, by = "initialPixels"]$N
      expect_true(all(frequ[attrib] <= floor(exactSizes[order(sams)][attrib])))
      expect_true(all(frequ[!attrib] == floor(exactSizes[order(sams)][!attrib])))
    } else {
      expect_true(all(out[, .N, by = "initialPixels"]$N == floor(exactSizes[order(sams)])))
    }
  }

  if (interactive()) print("testing exact maxSize, can't be achieved, allow jumping")
  exactSizes <- c(154, 111, 134) # too big for landscape, can't achieve it --
                                 #  will hit max numRetries, and will try jumping
  for (i in 1:20) {
    seed <- sample(1e6, 1)
    set.seed(seed)
    sams <- sample(innerCells, 3)
    out <- spread2(a, start = sams, spreadProb = 0.225, exactSize = exactSizes, asRaster = FALSE)
    expect_true(all(out[, .N, by = "initialPixels"]$N < exactSizes))
    expect_true(all(out$numRetries == 11)) # current max
  }

  if (interactive()) print("test circle = TRUE")
  for (i in 1:20) {
    message(i)
    seed <- sample(1e6, 1)
    set.seed(seed)
    sams <- sample(innerCells, length(sams))
    expect_error(spread2(a, start = sams, spreadProb = runif(1, 1.00000001, 1e4),
                         circle = TRUE, asRaster = FALSE, plot.it = TRUE))
    expect_error(spread2(a, start = sams, spreadProb = runif(1, -1e5, -0.00000001, 1e4),
                         circle = TRUE, asRaster = FALSE, plot.it = TRUE))
    out <- spread2(a, start = sams, spreadProb = 1, circle = TRUE, asRaster = FALSE)
    expect_true(is.numeric(out$distance))
    expect_true(NROW(out) == ncell(a))
  }

  # test circle
  sams <- sort(sample(innerCells, 3)) # sorted -- makes comparisons later easier
  out <- spread2(a, start = sams, spreadProb = 1, circle = TRUE, asRaster = FALSE,
                 returnDistances = TRUE)
  expect_true(NROW(out) == ncell(a))
  expect_true(all(out$state == "inactive"))
  expect_true(all(out$distance <= (sqrt(2) * ncol(a))))

  out <- spread2(a, start = sams, spreadProb = 1, circle = TRUE, allowOverlap = TRUE,
                 asRaster = FALSE, returnDistances = TRUE)
  expect_true(NROW(out) == ncell(a) * length(sams))
  expect_true(all(out$state == "inactive"))
  expect_true(all(out$distance <= (sqrt(2) * ncol(a))))

  setkey(out, initialPixels, distance)

  if (interactive()) {
    count <- 1
    for (ids in unique(out$initialPixels)) {
      dev(3 + count)
      count <- count + 1
      ras <- raster(a)
      ras[] <- 0
      ras[out[initialPixels == ids, pixels]] <- out[initialPixels == ids, distance]
      clearPlot()
      Plot(ras)
    }
  }

  if (interactive()) print("compare spread2 circle with cir circle")
  cirOut <- data.table(
    cir(a, allowOverlap = TRUE, loci = sams, minRadius = 0, maxRadius = 15,
        returnDistances = TRUE, simplify = TRUE)
  )
  if (interactive()) {
    for (ids in seq(unique(cirOut$id))) {
      dev(3 + ids)
      ras[cirOut[id == ids, indices]] <- cirOut[id == ids, dists]
      clearPlot()
      Plot(ras)
    }
  }
  cirOut$dists <- round(cirOut$dists, 4)
  out$distance <- round(out$distance, 4)
  setkey(cirOut, id, dists)
  #quickDT <- data.table(id = seq_along(sams), initialPixels = sams, key = "id")
  cirOut <- unique(cirOut)
  #cirOut <- quickDT[cirOut, on = ]
  setnames(cirOut, "id", "initialPixels")
  compare <- out[cirOut, on = c(initialPixels = "initialPixels", pixels = "indices")]
  expect_true(sum(abs(compare$dists - compare$distance)) %==% 0)

  if (interactive()) print("Scales with number of starts, not maxSize of raster")
  set.seed(21)
  b <- raster(extent(0, 33000, 0, 33000), res = 1)
  sams <- sample(ncell(b), 2)
  st1 <- system.time({
    out <- spread2(b, start = sams, spreadProb = 0.225, allowOverlap = TRUE, asRaster = FALSE)
  })
  expect_lt(st1[1], 1)

  if (interactive()) print("test neighProbs")
  maxSizes <- 14
  sp <- raster(a)
  spreadProbOptions <- 1:5
  sp[] <- sample(spreadProbOptions, ncell(sp), replace = TRUE)
  set.seed(2123)
  sams <- sample(innerCells, 2)
  set.seed(321)
  out <- spread2(a, spreadProb = 1, spreadProbRel = sp, start = sams,
                 neighProbs = c(0.7, 0.3), maxSize = maxSizes, asRaster = FALSE)
  expect_true(uniqueN(out) == maxSizes * length(sams))
  expect_true(NROW(out) == maxSizes * length(sams))

  if (interactive()) print("check variable lengths of neighProbs")
  set.seed(29937)
  sams <- sample(innerCells, 2)
  for (i in 1:8) {
    alwaysN <- rep(0, i)
    alwaysN[i] <- 1
    out <- spread2(
      a,
      spreadProb = 1,
      spreadProbRel = sp,
      iterations = 1,
      start = sams,
      neighProbs = alwaysN,
      asRaster = FALSE
    )
    expect_true(NROW(out) == (length(alwaysN) * 2 + length(sams)))
  }

  if (interactive()) {
    print(
      paste(
        "Test that when using neighProbs & a Raster of spreadProbs,",
        "the spreadProb raster is followed probabilistically.",
        "This test does only 1 iteration from 2 pixels that are",
        "not interacting with edges or each other"
      )
    )
  }
  sams <- sort(c(0:2 * 3 + 12) + rep(c(0, 30, 60, 90), 3))
  sams <- sams[sams < 90]
  set.seed(654)
  out <- list()
  for (i in 1:10) {
    out[[i]] <- spread2(a, spreadProbRel = sp, spreadProb = 1, iterations = 1,
                        start = sams, neighProbs = c(1), asRaster = FALSE)
  }
  out <- rbindlist(out)[state == "activeSource"]
  uniquePixels <- out[, list(uniquePix = unique(pixels)), by = "initialPixels"]
  avail <- table(sp[uniquePixels$uniquePix])
  actual <- unname(table(sp[out$pixels]))
  relProbs <- spreadProbOptions / sum(spreadProbOptions)
  aa <- rmultinom(1, size = 1e4, prob = relProbs)[, 1] * unname(avail)
  suppressWarnings(cht <- chisq.test(x = cbind(aa, actual)))

  if (as.numeric_version(paste0(R.version$major, ".", R.version$minor)) < "3.6.0") {
    expect_true(cht$p.value > 0.05)
  } else {
    expect_false(cht$p.value > 0.05) ## TODO: is this valid/correct test?
  }

  print("Scales with number of starts, not maxSize of raster")
  set.seed(21)
  b <- raster(extent(0, 10, 0, 10), res = 1)
  bProb <- gaussMap(b, speedup = 1)

  set.seed(1232)
  out <- spread2(spreadProb = 0.5, landscape = b, asRaster = FALSE,
                 start = ncell(b) / 2 - ncol(b) / 2, spreadProbRel = bProb,
                 returnFrom = TRUE, neighProbs = c(0.3, 0.7), exactSize = 30)

  set(out, NULL, "relProb", bProb[][out$pixels])
  if (interactive()) out

  if (interactive())
    print("check wide range of spreadProbs and that it makes a RasterLayer")
  set.seed(654)
  rasts <- list()
  for (i in 1:20) {
    rasts[[i]] <- spread2(a, spreadProb = stats::runif(1, 0, 1))
    expect_that(rasts[[i]], is_a("RasterLayer"))
  }
  if (interactive()) {
    names(rasts) <- paste0("ras", 1:20)
    clearPlot();Plot(rasts)
  }

  if (interactive())
    print("testing iterative calling of spread2")
  set.seed(299)
  sams <- sample(innerCells, 2)
  set.seed(299)
  out <- spread2(a, iterations = 1, start = sams, asRaster = FALSE)
  stillActive <- TRUE
  while (stillActive) {
    stillActive <- any(out$state == "activeSource")
    out <- spread2(a, iterations = 1, start = out, asRaster = FALSE)
  }
  set.seed(299)
  out2 <- spread2(a, start = sams, asRaster = FALSE)
  keyedCols <- c("initialPixels", "pixels")
  expect_equivalent(out2, out)

  if (interactive())
    print("testing iterative calling of spread2, but asRaster = TRUE")
  set.seed(299)
  sams <- sample(innerCells, 2)
  set.seed(299)
  out1 <- spread2(a, iterations = 1, start = sams, asRaster = TRUE)
  stillActive <- TRUE
  while (stillActive) {
    stillActive <- any(attr(out1, "pixel")$state == "activeSource")
    out1 <- spread2(a, iterations = 1, start = out1, asRaster = TRUE)
  }
  expect_true(identical(out, attr(out1, "pixel")))

  if (interactive())
    print("testing iterative with maxSize")
  set.seed(299)
  seed <- sample(1e6, 1)
  set.seed(seed)
  sams <- sample(innerCells, 2)
  exactSizes <- 5:6
  out <- spread2(a, start = sams, spreadProb = 0.225, iterations = 1,
                  exactSize = exactSizes, asRaster = FALSE)
  for (i in 1:20) {
    out <- spread2(a, start = out, spreadProb = 0.225, iterations = 1,
                   exactSize = exactSizes, asRaster = FALSE)
  }

  if (interactive())
    print("testing iterative with maxSize -- where needRetry occurs")
  set.seed(299)
  sams <- sample(innerCells, 2)
  exactSizes <- 60:61

  out <- spread2(a, start = sams, spreadProb = 0.225, iterations = 1,
                 exactSize = exactSizes, asRaster = FALSE)
  out2 <- spread2(a, start = sams, spreadProb = 0.225, iterations = 1,
                  exactSize = exactSizes, asRaster = FALSE)
  for (i in 1:25) {
    out <- spread2(a, start = out, spreadProb = 0.225, iterations = 1,
                   exactSize = exactSizes, asRaster = FALSE)
    attr(out2, "spreadState") <- NULL
    out2 <- spread2(a, start = out2, spreadProb = 0.225, iterations = 1,
                    exactSize = exactSizes, asRaster = FALSE)
  }
  expect_true(is.data.table(out))
  expect_true(is.data.table(out2))
  expect_true(all(attr(out2, "spreadState")$clusterDT$numRetries == 0))
  expect_true(all(attr(out, "spreadState")$clusterDT$numRetries > 10))

  # because loses info on how many retries, it will always be smaller
  expect_true(all(attr(out, "spreadState")$clusterDT$numRetries >
                    attr(out2, "spreadState")$clusterDT$numRetries))

  sams <- c(25, 75)
  set.seed(234)
  out <- spread2(a, start = sams, spreadProb = 0.225, iterations = 1,
                 exactSize = exactSizes, asRaster = FALSE)
  set.seed(234)
  out2 <- spread2(a, start = sams, spreadProb = 0.225, iterations = 1,
                  exactSize = exactSizes, asRaster = FALSE)
  for (i in 1:4) {
    # limit this so it doesn't get into retries, which will cause them to differ
    set.seed(234)
    out <- spread2(a, start = out, spreadProb = 0.225, iterations = 1,
                   exactSize = exactSizes, asRaster = FALSE)

    attr(out2, "spreadState") <- NULL
    set.seed(234)
    out2 <- spread2(a, start = out2, spreadProb = 0.225, iterations = 1,
                    exactSize = exactSizes, asRaster = FALSE)
  }

  ## they start to diverge if there is a jump that occurs, because the one without
  ## memory doesn't know how many retries it has had
  #expect_identical(data.table(out2), data.table(out)) ## TODO: fix this test

  for (i in 1:25) {
    set.seed(234)
    out <- spread2(a, start = out, spreadProb = 0.225, iterations = 1,
                   exactSize = exactSizes, asRaster = FALSE)

    attr(out2, "spreadState") <- NULL
    set.seed(234)
    out2 <- spread2(a, start = out2, spreadProb = 0.225, iterations = 1,
                    exactSize = exactSizes, asRaster = FALSE)
  }
  expect_false(identical(data.table(out2), data.table(out)))

  ##############################################################
  skip("benchmarking spread2")
  exactSizes <- 60:61
  microbenchmark(
    times = 10,
    a = {
      out <- spread2(a, start = sams, spreadProb = 0.225, iterations = 1,
                     exactSize = exactSizes, asRaster = FALSE)
      for (i in 1:25) {
        out <- spread2(a, start = out, spreadProb = 0.225, iterations = 1,
                       exactSize = exactSizes, asRaster = FALSE)
      }
    },
    b = {
      out2 <- spread2(a, start = sams, spreadProb = 0.225, iterations = 1,
                      exactSize = exactSizes, asRaster = FALSE)
      for (i in 1:25) {
        attr(out2, "spreadState") <- NULL
        out2 <- spread2(a, start = out2, spreadProb = 0.225, iterations = 1,
                        exactSize = exactSizes, asRaster = FALSE)
      }
    }
  )

  a <- raster(extent(0, 1000, 0, 1000), res = 1)
  set.seed(123)
  sams <- sample(innerCells, 30)
  set.seed(123)
  profvis::profvis({
    out <- spread2(a, start = sams, spreadProb = 0.235, asRaster = FALSE)
  })
  set.seed(123)
  profvis::profvis({
    out <- spread2(a, start = sams, spreadProb = 0.235, asRaster = FALSE, allowOverlap = TRUE)
  })

  set.seed(123)
  microbenchmark(times = 30, {
    out1 <- spread2(a, start = sams, spreadProb = 0.235, asRaster = FALSE)
  },
  b = {
    out2 <- spread(a, loci = sams, spreadProb = 0.235, id = TRUE)
  },
  c = {
    out2 <- spread(a, loci = sams, spreadProb = 0.235, id = TRUE, lowMemory = TRUE)
  })
  set.seed(123)
  profvis::profvis({
    out <- spread2(a, start = sams, spreadProb = 0.235, asRaster = FALSE, allowOverlap = TRUE)
  })

  ######## Benchmarking ##########
  iterativeFun <- function(a, skipChecks, n, sp) {
    sams <- sample(innerCells, n)
    out <- spread2(a, iterations = 1, start = sams, skipChecks = skipChecks,
                   asRaster = FALSE, spreadProb = sp)
    stillActive <- TRUE
    while (stillActive) {
      stillActive <- any(out$state == "activeSource")
      out <- spread2(a, spreadProb = sp, iterations = 1, start = out,
                     asRaster = FALSE, skipChecks = skipChecks)
    }
    out
  }

  nonIterativeFun <- function(a, skipChecks, n, sp) {
    sams <- sample(innerCells, n)
    out <- spread2(a, start = sams, asRaster = FALSE, skipChecks = skipChecks, spreadProb = sp)
    out
  }

  origSpread <- function(a, quick, n, sp) {
    sams <- sample(innerCells, n)
    out <- spread(a, spreadProb = sp, loci = sams, id = TRUE, returnIndices = TRUE,
                  quick = quick)
    out
  }

  origSpreadIterations <- function(a, quick, n, sp) {
    sams <- sample(innerCells, n)
    out <- spread(a, spreadProb = sp, iterations = 1, loci = sams, returnIndices = TRUE)
    stillActive <- TRUE
    while (stillActive) {
      stillActive <- any(out$active)
      out <- spread(a, spreadProb = sp, iterations = 1, spreadState = out,
                    returnIndices = TRUE, quick = quick)
    }
    out
  }

  n <- 2
  ras <- raster(extent(0, 160, 0, 160), res = 1)
  n <- 2000
  ras <- raster(extent(0, 6000, 0, 6000), res = 1)
  sp <- 0.225
  b <- raster(ras)
  b[] <- 1
  bb <- focal(b, matrix(1 / 9, nrow = 3, ncol = 3), fun = sum, pad = TRUE, padValue = 0)
  innerCells <- which(bb[] %==% 1)

  microbenchmark(
    times = 3,
    iterativeFun(ras, TRUE, n, sp),
    nonIterativeFun(ras, TRUE, n, sp),
    origSpread(ras, TRUE, n, sp),
    origSpreadIterations(ras, TRUE, n, sp)
  )
  #         iterativeFun(ras, TRUE, n, sp) 2.581752 2.667209 2.727886 2.752666 2.800953 2.849239
  #      nonIterativeFun(ras, TRUE, n, sp) 2.009914 2.268422 2.486540 2.526930 2.724854 2.922777
  #           origSpread(ras, TRUE, n, sp) 1.003594 1.043065 1.085511 1.082536 1.126470 1.170404
  # origSpreadIterations(ras, TRUE, n, sp) 7.267927 7.630435 7.908982 7.992943 8.229510 8.46607
#
  # without "skipChecks"
  microbenchmark(
    times = 100,
    iterativeFun(ras, FALSE, n, sp),
    nonIterativeFun(ras, FALSE, n, sp),
    origSpread(ras, FALSE, n, sp),
    origSpreadIterations(ras, FALSE, n, sp)
  )
  ## Unit: milliseconds
  ##                                    expr       min        lq      mean   median        uq       max neval
  ##         iterativeFun(ras, FALSE, n, sp)  3.096979 12.642477  73.02248 35.17520  91.10528  764.4073   300
  ##      nonIterativeFun(ras, FALSE, n, sp)  1.509484  6.555565  31.18444 14.91066  42.78317  158.5237   300
  ##           origSpread(ras, FALSE, n, sp)  5.154006  7.555631  14.87158 11.49005  17.50599  231.5487   300
  ## origSpreadIterations(ras, FALSE, n, sp) 10.754669 51.524368 141.48620 93.61996 169.10808 2110.2683   300
  ##
  profvis::profvis({
    set.seed(3451)
    for (i in 1:8)
      iterativeFun(ras, TRUE, n, sp = sp)
  })
  profvis::profvis({
    nonIterativeFun()
  })

  library(raster)
  library(data.table)
  library(fpCompare)
  n <- 2
  ras <- raster(extent(0, 160, 0, 160), res = 1)
  sp <- 0.225
  b <- raster(ras)
  b[] <- 1
  bb <-
    focal(
      b,
      matrix(1 / 9, nrow = 3, ncol = 3),
      fun = sum,
      pad = TRUE,
      padValue = 0
    )
  innerCells <- which(bb[] %==% 1)

  iterativeNeigh <- function(a, skipChecks, n, sp, exactSizes,
                           neighProbs, ...) {
    sams <- sample(innerCells, n)
    out <- spread2(a, iterations = 1, spreadProb = sp, start = sort(sams),
                   neighProbs = neighProbs, exactSize = exactSizes,
                   skipChecks = skipChecks, asRaster = FALSE, ...)
    stillActive <- TRUE
    while (stillActive) {
      out <- spread2(a, iterations = 1, spreadProb = sp, start = out,
                     neighProbs = neighProbs, exactSize = exactSizes,
                     skipChecks = skipChecks, asRaster = FALSE, ...)
      stillActive <- any(length(attr(out, "spreadState")$whActive) |
                           length(attr(out, "spreadState")$whNeedRetry))
    }
    out
  }

  origSpreadIterationsNeighs <- function(a, quick, n, sp, neighProbs, exactSize) {
    sams <- sample(innerCells, n)
    out <- spread(a, spreadProb = sp, neighProbs = neighProbs, iterations = 1,
                  loci = sams, exactSizes = TRUE, returnIndices = TRUE, maxSize = exactSize)
    stillActive <- TRUE
    while (stillActive) {
      stillActive <- any(out$active)
      out <- spread(a, spreadProb = sp, iterations = 1, neighProbs = neighProbs,
                    spreadState = out, maxSize = exactSize, exactSizes = TRUE,
                    returnIndices = TRUE, quick = quick)
    }
    out
  }

  sp <- randomPolygons(ras, numTypes = 35)
  sp[] <- (sp[] %% 5 + 1) / 10
  exactSizes <- c(123, 2240)
  neighProbs <- c(0.5, 0.3, 0.2)

  microbenchmark(
    times = 9,
    a = {
      iterativeNeigh(ras, TRUE, length(exactSizes), sp = sp, exactSize = exactSizes,
                     neighProbs = neighProbs)
    },
    b = {
      origSpreadIterationsNeighs(ras, TRUE, length(exactSizes), sp = sp,
                               exactSize = exactSizes, neighProbs = neighProbs)
    }
  )

  profvis::profvis({
    for (i in 1:5) {
      iterativeNeigh(ras, TRUE, length(exactSizes), sp = sp, exactSize = exactSizes,
                     neighProbs = neighProbs)
    }
  })

  profvis::profvis({
    for (i in 1:5) {
      origSpreadIterationsNeighs(ras, TRUE, length(exactSizes), sp = sp,
                                 exactSize = exactSizes, neighProbs = neighProbs)
    }
  })

  dev()
  iterativeNeigh(ras, TRUE, length(exactSizes), sp = sp,  exactSize = exactSizes, plot.it = TRUE)
  iterativeNeigh(ras, TRUE, length(exactSizes), sp = sp)#,  plot.it = TRUE)
  iterativeNeigh(ras, TRUE, length(exactSizes), sp = sp)#,  plot.it = TRUE)


  # compare original spread and spread2 -- seems pretty dead on
  nn <- 1000
  outNew <- out <- numeric(nn)
  for (i in 1:nn) {
    outNew[i] <- NROW(nonIterativeFun(ras, TRUE, n, sp))
    out[i] <- NROW(origSpread(ras, TRUE, n, sp))
  }

  out <- data.table(x = out)
  outNew <- data.table(x = outNew)
  mean(out$x)
  mean(outNew$x)
  sd(out$x)
  sd(outNew$x)

  n <- 5
  ras <- raster(extent(0, 1000, 0, 1000), res = 1)
  sp <- 0.295
  set.seed(123)
  microbenchmark(
    times = 100,
    nonIterativeFun(ras, TRUE, n, sp)
  )
})

test_that("spread2 tests -- asymmetry", {
  library(raster); on.exit(detach("package:raster"), add = TRUE)
  library(data.table); on.exit(detach("package:data.table"), add = TRUE)
  library(fpCompare); on.exit(detach("package:fpCompare"), add = TRUE)
  library(CircStats); on.exit(detach("package:CircStats"), add = TRUE)
  library(quickPlot); on.exit(detach("package:quickPlot"), add = TRUE)

  # inputs for x
  a <- raster(extent(0, 100, 0, 100), res = 1)
  b <- raster(a)
  b[] <- 1
  bb <- focal(b, matrix(1 / 9, nrow = 3, ncol = 3), fun = sum, pad = TRUE, padValue = 0)
  innerCells <- which(bb[] %==% 1)

  set.seed(123)
  sams <- sample(innerCells, 2)
  out <- spread2(a, start = sams, 0.215, asRaster = FALSE, asymmetry = 2,
                 asymmetryAngle = 90)
  for (i in 1:20) {
    expect_silent(
      out <- spread2(a, start = out, 0.215, asRaster = FALSE, asymmetry = 2,
                     asymmetryAngle = 90)
    )
  }

  a <- raster(extent(0, 1e2, 0, 1e2), res = 1)
  hab <- gaussMap(a, speedup = 1) # if raster is large (>1e6 pixels), use speedup>1
  names(hab) <- "hab"
  hab2 <- hab > 0
  maxRadius <- 25
  maxVal <- 50
  set.seed(53432)

  startCells <- as.integer(sample(1:ncell(hab), 1))

  n <- 16
  avgAngles <- numeric(n)
  lenAngles <- numeric(n)

  # function to calculate mean angle -- returns in degrees
  meanAngle <- function(angles) {
    deg(atan2(mean(sin(rad(angles))), mean(cos(rad(angles)))))
  }

  if (interactive()) {
    dev()
    clearPlot()
  }
  seed <- sample(1e6, 1)
  set.seed(seed)
  for (asymAng in (2:n)) {
    circs <- spread2(hab, spreadProb = 0.25, start = ncell(hab) / 2 - ncol(hab) / 2,
                    asymmetry = 40, asymmetryAngle = asymAng * 20, asRaster = FALSE)
    ci <- raster(hab)
    ci[] <- 0
    ci[circs$pixels] <- circs$initialPixels
    ciCentre <- raster(ci)
    ciCentre[] <- 0
    ciCentre[unique(circs$initialPixels)] <- 1
    newName <- paste0("ci", asymAng * 20)
    assign(newName, ci)

    where2 <- function(name, env = parent.frame()) {
      # simplified from pryr::where
      if (exists(name, env, inherits = FALSE)) env else where2(name, parent.env(env))
    }
    env <- where2(newName)
    if (interactive()) {
      objToPlot <- get(newName, envir = env)
      Plot(objToPlot, addTo = newName)
      Plot(ciCentre, cols = c("transparent", "black"), addTo = "objToPlot")
    }
    a <- cbind(id = circs$initialPixels, to = circs$pixels, xyFromCell(hab, circs$pixels))
    initialLociXY <- cbind(id = unique(circs$initialPixels),
                           xyFromCell(hab, unique(circs$initialPixels)))
    dirs <- directionFromEachPoint(from = initialLociXY, to = a)
    dirs[, "angles"] <- CircStats::deg(dirs[, "angles"])
    avgAngles[asymAng] <- tapply(dirs[, "angles"], dirs[, "id"], meanAngle) %% 360
    lenAngles[asymAng] <- tapply(dirs[, "angles"], dirs[, "id"], length)
  }

  whBig <- which(lenAngles > 50)
  pred <- (1:n)[whBig] * 20
  expect_true(abs(coef(lm(avgAngles[whBig] ~ pred))[[2]] - 1) < 0.1)

  # test that the events spread to the middle
  # Create a raster with one point at the centre
  ciCentre <- raster(hab)
  ciCentre <- setValues(ciCentre, 1)
  ciCentre[seq_len(ncell(ciCentre))[-(ncell(ciCentre) / 2 - ncol(ciCentre) / 2)]] <- NA_integer_
  # create a direction raster with all points leading to that point
  directionRas <- direction(ciCentre)
  directionRas[] <- deg(directionRas[])

  seed <- 4406
  set.seed(seed)
  sams <- ncol(directionRas) + 2
  circs <- spread2(hab, spreadProb = 0.265, start = sams, asymmetry = 300,
                   asymmetryAngle = directionRas, asRaster = TRUE)
  circs2 <- spread2(hab, spreadProb = 0.265, start = sams, asRaster = TRUE)
  if (interactive()) {
    Plot(circs, new = TRUE)
    ciCentrePlot <- ciCentre
    ciCentrePlot[ciCentrePlot == 2] <- NA
    ciCentrePlot[sams] <- 2
    Plot(ciCentrePlot, cols = c("transparent", "black", "red"), addTo = "circs")
    Plot(circs2, addTo = "circs", cols = "#1211AA33")
  }
  #test whether it stopped before hitting the whole map
  expect_true(sum(circs[], na.rm = TRUE) < ncell(circs))

  if (as.numeric_version(paste0(R.version$major, ".", R.version$minor)) < "3.6.0") {
    #test that it reached the centre, but not circs2 that did not have directionality
    expect_equal(circs[sams], circs[which(ciCentre[] == 1)]) ## TODO: restore this test
  }
  expect_true(is.na(circs2[ciCentre == 1]))
  expect_true(!is.na(circs2[sams]))

  # Here, test that the asymmetry version, with adjusted downward spreadProb is creating the
  #  same size events as the Non-asymmetry one. This is a weak test, really. It should
  sizes <- data.frame(a = numeric())
  set.seed(1234)
  for (i in 1:10) {
    sams <- ncell(hab) / 4 - ncol(hab) / 4 * 3
    circs <- spread2(hab, spreadProb = 0.18, start = sams,
                     asymmetry = 2, asymmetryAngle = 135, asRaster = TRUE)
    sizes <- rbind(sizes, cbind(a = attr(circs, "pixel")[, .N]))
    if (FALSE) {
      Plot(circs, new = TRUE)
      ciCentre[ciCentre == 2] <- NA
      ciCentre[sams] <- 2
      Plot(ciCentre, cols = c("black", "red"), addTo = "circs")
      Plot(circs2, addTo = "circs", cols = "#1211AA33")
    }
  }

  ttestOut <- t.test(sizes$a, mu = 994)
  expect_true(ttestOut$p.value > 0.05)

  skip("DEoptim within spread2")
  # This code is used to get the mean value for the t.test above
  n <- 100
  sizes <- integer(n)
  for (i in 1:n) {
    circs <- spread2(hab, spreadProb = 0.225,
                     start = ncell(hab) / 4 - ncol(hab) / 4 * 3,
                     asRaster = FALSE)
    sizes[i] <- circs[, .N]
  }
  goalSize <- mean(sizes)

  library(parallel)
  # only need 10 cores for 10 populations in DEoptim
  cl <- makeCluster(pmin(10, detectCores() - 2))
  parallel::clusterEvalQ(cl, {
    library(SpaDES.tools)
    library(raster)
    library(fpCompare)
  })

  objFn <- function(sp, n = 20, ras, goalSize) {
    sizes <- integer(n)
    for (i in 1:n) {
      circs <- spread2(ras, spreadProb = sp,
                       start = ncell(ras) / 4 - ncol(ras) / 4 * 3,
                       asymmetry = 2, asymmetryAngle = 135,
                       asRaster = FALSE)
      sizes[i] <- circs[, .N]
    }
    abs(mean(sizes) - goalSize)
  }
  aa <- DEoptim(objFn, lower = 0.2, upper = 0.23,
                control = DEoptim.control(
                  cluster = cl, NP = 10, VTR = 0.02,
                  initialpop = as.matrix(rnorm(10, 0.213, 0.001))
                ),
                ras = hab, goalSize = goalSize)

  # The value of spreadProb that will give the same expected event sizes to spreadProb = 0.225 is:
  sp <- aa$optim$bestmem
  circs <- spread2(ras, spreadProb = sp, start = ncell(ras) / 4 - ncol(ras) / 4 * 3,
                   asymmetry = 2, asymmetryAngle = 135, asRaster = FALSE)

  ####### Calibration curve
  skip("Calibration curves")
  n <- 500
  ras <- raster(extent(0, 1000, 0, 1000), res = 1)
  sp <- runif(n, 0.15, 0.25)
  sizes <- integer()
  for (i in 1:n) {
    circs <- spread2(ras, spreadProb = sp[i], start = ncell(ras) / 2 - ncol(ras) / 2,
                   asRaster = FALSE)
    sizes[i] <- NROW(circs)
    message(i)
  }
  dt1 <- data.table(sp, sizes)
  library(mgcv)
  aa <- gam(log10(dt1$sizes) ~ s(dt1$sp))
  aap <- predict(aa, se.fit = FALSE)
  plot(dt1$sp, log10(dt1$sizes), axes = FALSE, ylab = "Fire Size, ha", xlab = "Spread Probability")
  axis(2, 0:5, labels = 10 ^ (0:5))
  axis(1)
  aapOrd <- order(dt1$sp)
  lines(dt1$sp[aapOrd], aap[aapOrd], lwd = 2, col = "red")
  mtext(side = 3, paste("Resulting fire sizes, for given spread probabilities",
                        "Red line shows expected size", sep = "\n"))

  aa1 <- gam(dt1$sp ~ s(log10(dt1$sizes)))
  aap1 <- predict(aa1, se.fit = FALSE, type = "response")
  plot(log10(dt1$sizes), dt1$sp, axes = FALSE, xlab = "Fire Size, ha", ylab = "Spread Probability")
  axis(2)
  axis(1, 0:5, labels = 10 ^ (0:5))
  aap1Ord <- order(log10(dt1$sizes))
  lines(log10(dt1$sizes)[aap1Ord], aap1[aap1Ord], lwd = 2, col = "red")
  mtext(side = 3, paste("Resulting fire sizes, for given spread probabilities",
                        "Red line shows expected size", sep = "\n"))
})

test_that("spread2 returnFrom", {
  library(CircStats); on.exit(detach("package:CircStats"), add = TRUE)
  library(data.table); on.exit(detach("package:data.table"), add = TRUE)
  library(fpCompare); on.exit(detach("package:fpCompare"), add = TRUE)
  library(raster); on.exit(detach("package:raster"), add = TRUE)

  # inputs for x
  a <- raster(extent(0, 10, 0, 10), res = 1)
  b <- raster(a)
  b[] <- 1
  bb <- focal(b, matrix(1 / 9, nrow = 3, ncol = 3), fun = sum, pad = TRUE, padValue = 0)
  innerCells <- which(bb[] %==% 1)

  set.seed(123)
  for (i in 1:20) {
    sams <- sample(innerCells, 2)
    expect_silent(out <- spread2(a, start = sams, 0.215, asRaster = FALSE,
                                 returnFrom = TRUE))
    out <- spread2(a, start = sams, 0.215, asRaster = FALSE, returnFrom = TRUE)
    expect_true("from" %in% colnames(out))
    expect_true(sum(is.na(out$from)) == length(sams))
  }
})

test_that("spread2 tests", {
  library(raster)
  library(data.table)
  library(fpCompare)
  library(quickPlot)

  on.exit({
    detach("package:quickPlot")
    detach("package:fpCompare")
    detach("package:data.table")
    detach("package:raster")
  }, add = TRUE) # nolint

  # inputs for x
  a <- raster(extent(0, 100, 0, 100), res = 1)
  b <- raster(a)
  b[] <- 1
  bb <- focal(b, matrix(1 / 9, nrow = 3, ncol = 3), fun = sum, pad = TRUE, padValue = 0)
  innerCells <- which(bb[] %==% 1)
  sams <- sample(innerCells, 9)

  dev()
  expect_silent({
    out <- spread2(a, start = sams, 1, circle = TRUE, asymmetry = 4,
                   asymmetryAngle = 120, iterations = 10, asRaster = FALSE,
                   returnDistances = TRUE, allowOverlap = TRUE)
  })
  expect_true("effectiveDistance" %in% colnames(out))
  expect_true(all(out$state == "activeSource"))
  expect_true(all(out$distance[out$distance > 0] <= out$effectiveDistance[out$distance > 0]))
})

test_that("spread2 tests -- persistence", {
  library(raster)
  library(data.table)
  library(checkmate)
  library(bit)
  library(fastmatch)

  landscape <- raster::raster(nrows = 50, ncols = 50)
  landscape[] <- 1
  start <- 1:5

  ## test the effect of persistence as a single numeric value
  set.seed(5)
  noPersist <- spread2(landscape = landscape, start = start, asRaster = FALSE,
                       spreadProb = 0.23, persistProb = 0, iterations = 10, directions = 8L, plot.it = FALSE)
  wPersist <- spread2(landscape = landscape, start = start, asRaster = FALSE,
                      spreadProb = 0.23, persistProb = 0.8, iterations = 10, directions = 8L, plot.it = FALSE)

  expect_true(sum(noPersist$state == "activeSource") < sum(wPersist$state == "activeSource"))

  ## test the effect of persistence as a raster layer
  M <- matrix(0.8, nrow = 50, ncol = 50)
  M[upper.tri(M)] <- 0
  persistRas <- raster::raster(nrows = 50, ncols = 50)
  persistRas[] <- as.vector(M)

  ## first fire in high persistence area,
  ## second fire in low persistence area:
  start <- c(50, length(landscape)-49)

  set.seed(5)
  wRasPersist <- spread2(landscape = landscape, start = start,
                        spreadProb = 0.23, persistProb = persistRas, iterations = 10, directions = 8L,
                        asRaster = TRUE, plot.it = FALSE)

  expect_true(sum(wRasPersist[] == 1, na.rm = TRUE) > sum(wRasPersist[] == 2, na.rm = TRUE))

})


test_that("spread2 tests -- SpaDES.tools issue #22 NA in spreadProb", {
  library(raster)
  landscape <- raster::raster(nrows = 50, ncols = 50)
  landscape[] <- 1
  landscape[51:55] <- NA
  start <- 1:5
  spreadProb = landscape
  spreadProb[!is.na(landscape[])] <- runif(sum(!is.na(landscape[])))
  expect_silent(spread2(landscape = landscape, spreadProb = spreadProb,start = start,
          plot.it = FALSE))

})
