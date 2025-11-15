
# Create sample data mimicking CVI output
test_that("CVplot function works correctly", {
   
  Input1 <- data.frame(
  J1 = c(5, 5, 6, 6, 6, 6, 6, 6, 5, 6, 5, 6, 3, 4, 4),
  J2 = c(5, 2, 6, 6, 6, 6, 6, 5, 4, 6, 5, 6, 4, 4, 3),
  J3 = c(5, 5, 6, 6, 6, 5, 6, 5, 4, 5, 5, 6, 5, 3, 5),
  J4 = c(5, 5, 6, 6, 6, 6, 6, 6, 5, 6, 3, 6, 5, 3, 5),
  J5 = c(5, 5, 6, 6, 6, 6, 6, 6, 6, 5, 5, 6, 3, 4, 5))
  Input1.CVIoutput <- CVI(data = Input1, cut = 4, conf.level = .90)

  # Test basic functionality
  plot <- CVplot(data = Input1.CVIoutput, item.col = "Item", point.coeficient = "CVI",
  lwr.ci = "lwr.ci", up.ci = "upr.ci")

  # Test output
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot, c("gg", "ggplot"))
})


















