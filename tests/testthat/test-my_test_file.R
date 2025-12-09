

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})




test_that("tau_1D calculates shear stress correctly", {
  R <- 3.1
  S <- 0.003
  t <- 9.81 * 1000 * R * S
  expect_equal(tau_1D(R,S), t, tolerance = 0.001)
})

test_that("tau_1D calculates shear stress correctly with custom density and gravity", {
  R <- 3.1
  S <- 0.003
  gc <- 4
  rc <- 1200
  t <- 4 * 1200 * R * S
  expect_equal(tau_1D(R,
                      S,
                      g = gc,
                      rho = rc
                      ),
               t,
               tolerance = 0.001
               )
})

test_that("u_ferg calculates velocity correctly", {
  R <- 3.27
  S <- 0.00066
  d84 <- 0.02
  u <- 2.21  #taken from flood hydraulics tab of BGC scour
  expect_equal(u_ferg(R,S,d84),
               u,
               tolerance = 0.001
               )
})

test_that("u_ferg calculates velocity correctly with custom a1 and a2 values", {
  R <- 3.60
  S <- 0.00066
  d84 <- 0.02
  a1c = 5.5
  a2c = 2.20
  u <- 1.99  #taken from flood hydraulics tab of BGC scour (add path to..)
  expect_equal(u_ferg(R,
                     S,
                     d84,
                     a1 = a1c,
                     a2 = a2c
                     ),
               u,
               tolerance = 0.005
               )
})

test_that("t_star calculates dimesionless shear stress correctly", {
  tau <- 150
  d <- 0.12
  ts <- 0.0772
  #ts based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(t_star(tau,d),
               ts,
               tolerance = 0.001
               )
})

test_that("t_star calculates dimesionless shear stress correctly for custom densities and gravity", {
  tau <- 150
  d <- 0.12
  gc <- 4
  rc <- 1200
  rsc <- 2650
  ts <- 0.216
  #ts based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(t_star(tau,
                     d,
                     g = gc,
                     rho = rc,
                     rho_s = rsc
                     ),
               ts,
               tolerance = 0.005
               )
})

test_that("qb_conversion calculates sediment flux correctly", {
  q_star <- 0.251
  d <- 0.06
  qb <- 0.0148
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(qb_conversion(q_star,d),
               qb,
               tolerance = 0.005
               )
})

test_that("qb_conversion calculates sediment flux correctly with custom values", {
  q_star <- 1.095
  d <- 0.06
  gc <- 6
  Gsc <- 1.3
  qb <- 0.0449
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(qb_conversion(q_star,
                            d,
                            g = gc,
                            Gs = Gsc
                            ),
               qb,
               tolerance = 0.005
               )
})

test_that("t_crit_soul calculates tc correctly", {
  d <- 0.003
  tc <- 0.0462
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(t_crit_soul(d),
               tc,
               tolerance = 0.005
               )
})

test_that("t_crit_soul calculates tc correctly with custom inputs", {
  d <- 0.003
  gc <- 8
  rc <- 1000
  rsc <- 2550
  Gsc <- 1.55
  kvc <- 1.8e-6
  tc <- 0.0387
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(t_crit_soul(d,
                          g = gc,
                          rho = rc,
                          rho_s = rsc,
                          kv = kvc
                          ),
               tc,
               tolerance = 0.005
               )
})

test_that("t_crit_vr calculates tc correctly for 0.1 mm", {
  d <- 0.0001
  tc <- 0.0949
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(t_crit_vr(d),
               tc,
               tolerance = 0.005
  )
})

test_that("t_crit_vr calculates tc correctly for 0.3 mm", {
  d <- 0.0003
  tc <- 0.0383
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(t_crit_vr(d),
               tc,
               tolerance = 0.005
  )
})

test_that("t_crit_vr calculates tc correctly for 0.6 mm", {
  d <- 0.0006
  tc <- 0.0305
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(t_crit_vr(d),
               tc,
               tolerance = 0.005
  )
})
test_that("t_crit_vr calculates tc correctly for 3 mm", {
  d <- 0.003
  tc <- 0.0456
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(signif(t_crit_vr(d),
                      digits = 3
                      ),
               tc
               )
})

test_that("t_crit_vr calculates tc correctly for 30 mm", {
  d <- 0.03
  tc <- 0.056
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(t_crit_vr(d),
               tc,
               tolerance = 0.005
               )
})

test_that("t_crit_vr calculates tc correctly with custom inputs", {
  d <- 0.003
  gc <- 8
  rc <- 1000
  rsc <- 2550
  kvc <- 1.8e-6
  tc <- 0.0397
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(t_crit_vr(d,
                        g = gc,
                        rho = rc,
                        rho_s = rsc,
                        kv = kvc
                        ),
               tc,
               tolerance = 0.005
               )
})

test_that("d_star calculates d* correctly with custom inputs", {
  d <- 0.003
  gc <- 8
  Gsc <- 1.55
  kvc <- 1.8e-6
  ds <- 46.9
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(signif(d_star(d,
                             g = gc,
                             Gs = Gsc,
                             kv = kvc
                             ),
                      digits = 3
                      ),
               ds
               )
})

test_that("d_star calculates d* correctly", {
  d <- 0.003
  ds <- 75.9
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(d_star(d),
               ds,
               tolerance = 0.005
               )
})

test_that("mpm calculates q* correctly", {
  d <- 0.05
  R <- 1.7
  S <- 0.003
  qs <- 0.0144
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(mpm(R,
                  S,
                  d
                  ),
               qs,
               tolerance = 0.005
               )
})

test_that("mpm returns 0 for t < tc", {
  d <- 0.05
  R <- 1.7
  S <- 0.0003
  qs <- 0.0144
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_true(mpm(R, S, d) == 0)
})

test_that("mpm calculates q* correctly with custom densities", {
  d <- 0.05
  R <- 1.7
  S <- 0.003
  rc <- 1100
  rsc <- 2550
  qs <- 0.0424
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(mpm(R,
                  S,
                  d,
                  rho = rc,
                  rho_s = rsc
                  ),
               qs,
               tolerance = 0.005
               )
})

test_that("mpm calculates q* correctly with soulsby tc", {
  d <- 0.05
  R <- 1.7
  S <- 0.003
  tcc <- t_crit_soul(d)
  qs <- 0.00431
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(mpm(R,
                  S,
                  d,
                  t_crit = tcc
                  ),
               qs,
               tolerance = 0.005
               )
})

test_that("wp calculates q* correctly", {
  d <- 0.05
  R <- 1.7
  S <- 0.003
  qs <- signif(0.01443/2, digits = 3)
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(wp(R,S,d),
               qs,
               tolerance = 0.005
               )
})

test_that("wp returns 0 for t < tc", {
  d <- 0.05
  R <- 1.7
  S <- 0.0003
  qs <- 0.0144
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_true(wp(R, S, d) == 0)
})

test_that("wp calculates q* correctly with soulsby tc", {
  d <- 0.05
  R <- 1.7
  S <- 0.003
  tcc <- t_crit_soul(d)
  qs <- signif(0.00431 /2, digits = 3)
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(wp(R,
                 S,
                 d,
                 t_crit = tcc
                 ),
               qs,
               tolerance = 0.005
               )
})

test_that("vr calculates q* correctly", {
  d <- 0.001
  R <- 1.4
  S <- 0.0003
  qs <- 1.08
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(vr(R,
                 S,
                 d
                 ),
               qs,
               tolerance = 0.005
               )
})

test_that("vr calculates q* correctly with soulsby tc", {
  d <- 0.001
  R <- 1.4
  S <- 0.0003
  tcc <- t_crit_soul(d)
  qs <- 1.23
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(vr(R,
                 S,
                 d,
                 t_crit = tcc
                 ),
               qs,
               tolerance = 0.005
               )
})

test_that("vr calculates q* correctly with custom variables", {
  d <- 0.001
  R <- 1.4
  S <- 0.0003
  gc <- 7
  rc <- 1050
  rsc <- 2600
  kvc <- 1.3e-6
  tcc <- t_crit_soul(d,
                     rc,
                     rsc,
                     gc,
                     kvc
                     )
  qs <- 1.99
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(vr(R,
                 S,
                 d,
                 t_crit = tcc,
                 g = gc,
                 rho = rc,
                 rho_s = rsc,
                 kv = kvc
                 ),
               qs,
               tolerance = 0.005
               )
})

test_that("yln calculates q* correctly", {
  d <- 0.001
  R <- 1.4
  S <- 0.0003
  qs <- 1.05
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(yln(R, S, d),
               qs,
               tolerance = 0.005
               )
})
test_that("yln calculates q* correctly with Van Rijn tc", {
  d <- 0.001
  R <- 1.4
  S <- 0.0003
  tcc <- t_crit_vr(d)
  qs <- 0.967
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(yln(R,
                  S,
                  d,
                  t_crit = tcc
                  ),
               qs,
               tolerance = 0.001
               )
})
test_that("yln calculates q* correctly with custom variables", {
  d <- 0.001
  R <- 1.4
  S <- 0.0003
  gc <- 7
  rc <- 1050
  rsc <- 2600
  kvc <- 1.3e-6
  tcc <- t_crit_vr(d,
                   rc,
                   rsc,
                   gc,
                   kvc
                   )
  qs <- 1.44
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(yln(R,
                    S,
                    d,
                    t_crit = tcc,
                    g = gc,
                    rho = rc,
                    rho_s = rsc,
                    kv = kvc
                    ),
               qs,
               tolerance = 0.001
               )
})

test_that("eb calculates q* correctly for t* < 0.19", {
  d <- 0.03
  R <- 1.9
  S <- 0.003
  qs <- 0.0588
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(eb(R, S, d), qs, tolerance = 0.001)
})

test_that("eb calculates q* correctly for t* > 0.19", {
  d <- 0.03
  R <- 1.9
  S <- 0.006
  qs <- 0.399
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(eb(R, S, d), qs, tolerance = 0.001)
})

test_that("eb calculates q* correctly with custom inputs", {
  d <- 0.03
  R <- 1.9
  S <- 0.003
  gc <- 9.0
  rc <- 1100
  rsc <- 2700
  kvc <- 1.2e-6
  qs <- 0.0880
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_equal(eb(R,
                  S,
                  d,
                  g = gc,
                  rho = rc,
                  rho_s = rsc,
                  kv = kvc
                  ),
               qs,
               tolerance = 0.001
               )
})
test_that("ec calculates q* smaller than mpm", {
  d <- 0.01
  R <- 1.7
  S <- 0.001
  #ans based on online calculator at https://personalpages.manchester.ac.uk/staff/david.d.apsley/hydraulics/bedload.htm
  expect_lt(ec(R, S, d), mpm(R, S, d))
})

test_that("ec calculates q* similar to wp at high transport values", {
  d <- 0.01
  R <- 1.7
  S <- 0.01
  a <- log10(ec(R, S, d))
  b <- log10(wp(R,S, d))
  expect_lt(abs(a - b), 0.5)
})

test_that("ec calculates q* similar to eb at moderate transport values", {
  d <- 0.01
  R <- 0.7
  S <- 0.01
  a <- log10(ec(R, S, d))
  b <- log10(eb(R,S, d))
  expect_lt(abs(a - b), 0.5)
})

test_that("ec calculates q* similar to wp at low transport values", {
  d <- 0.01
  R <- 1.7
  S <- 0.001
  a <- log10(ec(R, S, d))
  b <- log10(wp(R,S, d))
  expect_lt(abs(a - b), 0.5)
})

test_that("ec calculates q* of 0 for very small stresses", {
  d <- 0.01
  R <- 1.7
  S <- 0.0001
  expect_true(ec(R,S,d) == 0)
})

test_that("Increases in fluid density have similar effects on q* for ec and wp for moderate transport rates", {
  d <- 0.01
  R <- 1.7
  S <- 0.001
  r1 <- 1000
  r2 <- 1300
  a <- ec(R, S, d, rho = r2) / ec(R, S, d, rho = r1)
  b <- wp(R, S, d, rho = r2) / wp(R, S, d, rho = r1)
  ndiff <- abs(a - b) / a
  expect_lt(ndiff, 0.25)
})

test_that("Increases in fluid density have similar effects on q* for ec and wp for high transport rates", {
  d <- 0.01
  R <- 1.7
  S <- 0.01
  r1 <- 1000
  r2 <- 1300
  a <- ec(R, S, d, rho = r2) / ec(R, S, d, rho = r1)
  b <- wp(R, S, d, rho = r2) / wp(R, S, d, rho = r1)
  ndiff <- abs(a - b) / a
  expect_lt(ndiff, 0.25)
})

test_that("Increases in sediment density have similar effects on q* for ec and wp for high transport rates", {
  d <- 0.01
  R <- 1.7
  S <- 0.01
  r1 <- 2600
  r2 <- 3000
  a <- ec(R, S, d, rho_s = r2) / ec(R, S, d, rho_s = r1)
  b <- wp(R, S, d, rho_s = r2) / wp(R, S, d, rho_s = r1)
  ndiff <- abs(a - b) / a
  expect_lt(ndiff, 0.25)
})

test_that("Increases in sediment density have similar effects on q* for ec and wp for moderate transport rates", {
  d <- 0.01
  R <- 1.7
  S <- 0.001
  r1 <- 2600
  r2 <- 3000
  a <- ec(R, S, d, rho_s = r2) / ec(R, S, d, rho_s = r1)
  b <- wp(R, S, d, rho_s = r2) / wp(R, S, d, rho_s = r1)
  ndiff <- abs(a - b) / a
  expect_lt(ndiff, 0.25)
})

test_that("q* values predicted using ec with and without specifying U should be the same if d84 <- 2 * d", {
  d <- 0.01
  R <- 1.7
  S <- 0.001
  Uf <- u_ferg(R,S, 2*d)
  a <- ec(R, S, d)
  b <- ec(R, S, d, U = Uf)
  expect_equal(a, b)
})

test_that("Increasing d84 in ec will decrease the transport rate, q* by reducing stream power", {
  d <- 0.01
  R <- 1.7
  S <- 0.001
  a <- ec(R, S, d, d84 = 2*d)
  b <- ec(R, S, d, d84 = 3*d)
  expect_gt(a, b)
})

test_that("Increasing d in ec will decrease the transport rate, q* by increasing critical stress", {
  d <- 0.01
  R <- 1.7
  S <- 0.001
  a <- ec(R, S, d, d84 = 2*d)
  b <- ec(R, S, 1.2*d, d84 = 2*d)
  expect_gt(a, b)
})

test_that("rk and ec produces similar qb for low stresss", {
  d <- 0.05
  d84 <- 2*d
  R <- 0.9
  S <- 0.005
  a <- qb_conversion(ec(R, S, d, d84),
                     d)
  b <- qb_conversion(rk(R, S, d, d84),
                     d84)
  dif <- abs(a - b) / a
  expect_equal(dif, 0, tolerance = 0.25)
})

test_that("rk and ec produces similar qb for high stresss", {
  d <- 0.05
  d84 <- 2*d
  R <- 0.9
  S <- 0.05
  a <- qb_conversion(ec(R, S, d, d84),
                     d)
  b <- qb_conversion(rk(R, S, d, d84),
                     d84)
  dif <- abs(a - b) / a
  expect_equal(dif, 0, tolerance = 0.25)
})


test_that("rk and ec produces similar qb for high stresss", {
  d <- 0.1
  d84 <- 2*d
  R <- 0.9
  S <- 0.05
  a <- qb_conversion(ec(R, S, d, d84),
                     d)
  b <- qb_conversion(rk(R, S, d, d84),
                     d84)
  dif <- abs(a - b) / a
  expect_equal(dif, 0, tolerance = 0.75)
})

test_that("rk and ec produces similar qb for high stresss", {
  d <- 0.05
  d84 <- 2*d
  R <- 0.6
  S <- 0.05
  a <- qb_conversion(ec(R, S, d, d84),
                     d)
  b <- qb_conversion(rk(R, S, d, d84),
                     d84)
  dif <- abs(a - b) / a
  expect_equal(dif, 0, tolerance = 0.5)
})

test_that("rk and ec produces similar qb for high stresss", {
  d <- 0.05
  d84 <- 2*d
  R <- 0.6
  S <- 0.005
  a <- qb_conversion(ec(R, S, d, d84),
                     d)
  b <- qb_conversion(rk(R, S, d, d84),
                     d84)
  dif <- abs(a - b) / a
  expect_equal(dif, 0, tolerance = 0.25)
})
