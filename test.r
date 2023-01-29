# Import data manually from Python output
a_pos <- 45593
a_tot <- 81885000
b_pos <- 8238
b_tot <- 6944237

#Chisquare
results <- c(a_pos, a_tot - a_pos, b_pos, b_tot - b_pos)
m <- matrix(results, nrow=2)
chisq.test(m, correct = F)

#Binomial
prop.test(c(b_pos, a_pos), c(b_tot, a_tot))

rr <- (b_pos / b_tot) / (a_pos / a_tot) # Relative risk 
se <- ((1 / b_pos) - (1 / b_tot) + (1 / a_pos) - (1 / a_tot)) # SE

print(rr) # Ratios
sprintf("RR: %f | CI 95%%: %f - %f", rr, rr - 1.96 * se, rr + 1.96 * se)
print(1.96 * se)
