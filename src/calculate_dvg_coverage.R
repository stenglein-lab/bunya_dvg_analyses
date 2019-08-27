library (stringr)
library (dplyr)
library (ggplot2)
library (gridExtra)

# this script calculates the relative abundance of particular DVGs (relative to median total coverage for a segment)
# and creates a couple plots of these abundances
#
# Mark Stenglein May 14, 2019

# import and tidy the data
if (!exists("depth_df")) source("./import_and_tidy_dvg_data.R")

# calculated median depth of (total) coverage for each segment and store it in a new df
median_depths <- depth_df %>% group_by(segment) %>% summarize(median_depth = median(depth))

# filter to keep deletion DVGs supported by > 3 reads and with a deletion > 1 bp
# these should've already been filtered by DI-tector command line parameters
del_df <- filter(del_df, counts >= 3, delta > 1)

# merge deletion DVG df with median depth df
del_df <- left_join(del_df, median_depths, by="segment")

# create new column bp_abundance with calculated relative breakpoint abundances: counts
# for that breakpoint divided by median coverage depth for that segment
del_df <- del_df %>% mutate(bp_abundance = counts / median_depth)

# sort by decreasing abundance
# this just outputs to console, doesn't actually sort df
del_df %>% arrange(-bp_abundance)

# calculate mean abundances for each type of segment (S, M, L), first removing NA values which mess up mean calculation
# this just outputs to console
del_df %>% filter(!is.na(bp_abundance)) %>% group_by(seg) %>% summarize(mean = mean(bp_abundance))

# how many deletion DVGs by each segment type (S, M, L)
# this just outputs to console
del_df %>% group_by(seg) %>% summarize(n=n())

# how many segments had any deletion DVGs?
num_dvgs_df <- del_df %>% group_by(segment) %>% summarize(num_dvgs = n())
nrow(num_dvgs_df)

# how many segments total?
num_segs_df <- depth_df %>% group_by(segment) %>% summarize(n=n())
nrow(num_segs_df)

# do some math
85 / 306
1108/85

# how many had relative abundances > 1% or 10%?
gt_1pct <- del_df %>% filter(bp_abundance > 0.01)
gt_10pct <- del_df %>% filter(bp_abundance > 0.1)


#
# make figure plotting relative DVG abundance vs. segment length
# color by segment type (S, M, L)
#
# note that we can calculate segment length from the DI-tector output, since
# DI-tector tells us the DVG length plus the length of the deletion
#

color_blind_friendly_color_trio <- c("#0072B2", "#009E73", "#F0E442")

# make the plot, save to an object
p <- ggplot(data=del_df, aes(x=(delta+length), y=bp_abundance)) +
  geom_point(aes(fill=seg), shape=21, color = "black", size = 2, stroke = 0.2, alpha=0.5) +
  xlab("Segment length (nt)") + 
  ylab("Normalized DVG abundance\n(supporting split reads / median depth)") +
  scale_y_log10() +
  scale_fill_manual(values=color_blind_friendly_color_trio)

# this will show the plot in the plot viewer in R studio
p

# save to a PDF
ggsave("abundance.pdf", p, height = 3 , width = 4 )


# this plot slightly different
# now grouped by S, M, L
# jitter plot
# and box plot overlay

p2 <- ggplot(data=del_df, aes(x=(seg), y=bp_abundance)) +
  geom_jitter(aes(fill=seg), shape=21, color = "black", size = 2, stroke = 0.2, alpha=0.75) +
  # geom_violin(aes(fill=seg), alpha=0.25) +
  geom_boxplot(aes(fill=seg), alpha=0.25) +
  xlab("Segment") +
  ylab("Normalized DVG abundance (split reads / median depth)") +
  scale_y_log10() +
  scale_fill_manual(values=color_blind_friendly_color_trio)

p2

