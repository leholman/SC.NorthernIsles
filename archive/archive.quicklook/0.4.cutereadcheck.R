### 
dat2 <- read_excel("metadata.xlsx")

barplot(dat2$)

dat2$rawreads <- unname(as.numeric(gsub("[^0-9.]", "", dat2$`raw reads`)) * c(N=1, K=1e3, M=1e6)[sub("^$", "N", toupper(gsub("[^A-Za-z]", "", dat2$`raw reads`)))])
dat2 <- dat2[order(dat2$CoreID,dat2$Est.Cal.YrBP), ]
unique(dat2$CoreID)
barplot(dat2$rawreads, col="grey80", border=NA, ylab="reads")
barplot(dat2$`premapping reads`, col="steelblue", border=NA, axes=FALSE, add=TRUE)


CORE_COLS <- c(GC017 ="#D55E00", 
               PC012 = "#0072B2", 
               PC014 = "#009E73")
cols <- CORE_COLS[dat2$CoreID]
cols[is.na(cols)] <- "grey"

big   <- dat2$rawreads / 1e6
small <- dat2$`premapping reads` / 1e6
top   <- max(big, na.rm = TRUE)

bp <- barplot(big, col = "grey80", border = NA, axes = FALSE,
              ylab = "Reads (millions)", ylim = c(-top * 0.10, top))
barplot(small, col = "steelblue", border = NA, axes = FALSE, add = TRUE)
axis(2, at = pretty(c(0, top)),las=1)

w <- diff(bp)[1] / 2
rect(bp - w, -top * 0.05, bp + w, -top * 0.02, col = cols, border = NA, xpd = NA)

