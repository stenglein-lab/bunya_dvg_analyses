library (stringr)
library (dplyr)
library (ggplot2)
library (ggExtra)
library (scales)
library (ggpubr)
library (gridExtra)
# library (cowplot)

# This script creates summary plots for all the deletion DVGs
#
# Mark Stenglein, May 3, 2019

# import and tidy the data
if (!exists("del_df")) source("./import_and_tidy_dvg_data.R")

# scatter plot of breakpoint positions
# (not used in paper figure)
sp <- ggplot() + geom_point(data=filter(del_df, delta>1), aes(x=bp, y=dist_3p_end, color=seg), alpha=0.5) + 
  xlab ("breakpoint distance from 5' end") + 
  ylab ("breakpoint distance from 3' end") +
  scale_colour_manual(values = c("navyblue", "red", "springgreen1")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_line(color="slategray1", linetype=5, size=0.2)) + 
  coord_fixed()
sp

# add histograms to the margins of the scatter plot
ggMarginal(sp, type = "histogram", binwidth=50, groupColour=T, groupFill=T)

# calculate some summary stats
# by S, M, or L
del_df %>% group_by (seg) %>% summarize(mean_bp = mean(bp), median_bp = median(bp), mean_3p_dist = mean(dist_3p_end), median_3p_dist = median(dist_3p_end))
# total
del_df %>% summarize(mean_bp = mean(bp), median_bp = median(bp), mean_3p_dist = mean(dist_3p_end), median_3p_dist = median(dist_3p_end))

# order the deletion DVGs by their bp positions
# TODO: Make this a function

# plot cartoons of all DVGs
# L segment first
# there are 800 some L segment DVGs

# figure out how tall graph should be and how thick lines should be to fit
dvgs_to_show <- 900
  
# the height of the plot in mm (1 inch ~ 25 mm)
total_plot_height <- 100

# estimate ~75% of the plot area is actual data (bars), rest title, axes, etc
data_plot_height <- total_plot_height * 0.75  

# this is the width of the line in mm
line_width <- (data_plot_height / dvgs_to_show)

# prepare L seg DVGs for plotting
cp_df_L <- del_df %>% filter(delta>1, counts >= 3, seg=="L") %>%        # keep only L segment DVGs with >1 bp deleted and >3 supporting reads 
                      filter(length+delta < 8000) %>%                   # keep only orthobunyavirus L segments, which have lengths < 8000
                      arrange(-counts) %>% filter(row_number() <= dvgs_to_show)    # keep only the N-most abundant DVGs
cp_df_L <- cp_df_L %>% arrange(delta) %>% mutate(plot_order=row_number())          # sort by decreasing deletion length.  plot_order = new column with order for plotting purposes (will by y variable).

# cartoon plot of DVGs
#
# the 3 geom_segment geometries make the orange boxes and the grey connecting line
#
cp_L <- ggplot(data=cp_df_L) +
  geom_segment(aes(x=1,  xend=bp, y=plot_order, yend=plot_order), size=line_width, color="darkorange") +
  geom_segment(aes(x=ri, xend=(length+delta), y=plot_order, yend=plot_order),  size=line_width, color="darkorange") +
  # geom_segment(aes(x=bp, xend=ri, y=plot_order, yend=plot_order, color=counts), size=line_width)  +
  geom_segment(aes(x=bp, xend=ri, y=plot_order, yend=plot_order), size=line_width, color="grey90")  +
  scale_color_gradient(low='grey90',high='blue', limits=c(0,20), oob=squish) +
  xlab ("genome position (nt)") + 
  ylab ("200 best supported DVGs") +
  scale_y_continuous(limits=c(1,dvgs_to_show))+
  theme(text = element_text(size=9, family="Helvetica"),
        axis.text.x = element_text(size=8, family="Helvetica"),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank())

# show in R studio viewer
cp_L

# M/S segments now

# figure out how tall graph should be and how thick lines should be to fit
# There are just less than 110 each M/S DVGs
dvgs_to_show <- 110
  
# the height of the plot in mm (1 inch ~ 25 mm)
# these will each be ~1/2 the height of the L plot
total_plot_height <- 50

# estimate ~75% of the plot area is actual data (bars), rest title, axes, etc
data_plot_height <- total_plot_height * 0.75  

# this is the width of the line in mm
line_width <- (data_plot_height / dvgs_to_show)

# reorder M, S DVGs
cp_df_M <- del_df %>% filter(delta>1, counts >= 3, seg=="M") %>% arrange(-counts) %>% filter(row_number() <= dvgs_to_show) 
cp_df_M <- cp_df_M %>% arrange(delta) %>% mutate(plot_order=row_number())

cp_df_S <- del_df %>% filter(delta>1, counts >= 3, seg=="S") %>% arrange(-counts) %>% filter(row_number() <= dvgs_to_show) 
cp_df_S <- cp_df_S %>% arrange(delta) %>% mutate(plot_order=row_number())

# cartoon plot of DVGs
cp_M <- ggplot(data=cp_df_M) +
  geom_segment(aes(x=1,  xend=bp, y=plot_order, yend=plot_order), size=line_width, color="darkorange") +
  geom_segment(aes(x=ri, xend=(length+delta), y=plot_order, yend=plot_order),  size=line_width, color="darkorange") +
  # geom_segment(aes(x=bp, xend=ri, y=plot_order, yend=plot_order, color=counts), size=line_width)  +
  geom_segment(aes(x=bp, xend=ri, y=plot_order, yend=plot_order), size=line_width, color="grey90")  +
  scale_color_gradient(low='grey90',high='blue', limits=c(0,20), oob=squish) +
  xlab ("genome position (nt)") + 
  ylab ("200 best supported DVGs") +
  scale_y_continuous(limits=c(1,dvgs_to_show))+
  theme(text = element_text(size=9, family="Helvetica"),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size=8, family="Helvetica"),
        axis.text.y = element_blank())
   
# cartoon plot of DVGs
cp_S <- ggplot(data=cp_df_S) +
  geom_segment(aes(x=1,  xend=bp, y=plot_order, yend=plot_order), size=line_width, color="darkorange") +
  geom_segment(aes(x=ri, xend=(length+delta), y=plot_order, yend=plot_order),  size=line_width, color="darkorange") +
  # geom_segment(aes(x=bp, xend=ri, y=plot_order, yend=plot_order, color=counts), size=line_width)  +
  geom_segment(aes(x=bp, xend=ri, y=plot_order, yend=plot_order), size=line_width, color="grey90")  +
  scale_color_gradient(low='grey90',high='blue', limits=c(0,20), oob=squish) +
  xlab ("genome position (nt)") + 
  scale_y_continuous(limits=c(1,dvgs_to_show))+
  # ylab ("200 best supported DVGs") +
  theme(text = element_text(size=9, family="Helvetica"),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size=8, family="Helvetica"),
        axis.text.y = element_blank())

# arrange the 3 plots in a grid: L on left, M and S on right
all_3_cp <- grid.arrange(cp_L, cp_M, cp_S, ncol=2, nrow=2, layout_matrix = rbind(c(1,2),c(1,3)))

# save them
ggsave("dvg_cartoons.pdf", plot=all_3_cp, height=100, width=75, units="mm")


# make histograms of breakpoints positions
# this is actually 2 superimposed histograms: one for the deletion bps (initiating breakpoint) and one for the deletion ris (re-initating breakpoints)
# (not used in paper figure)
histo_L <- ggplot(data = del_df %>% filter(delta>1, seg=="L") %>% filter(length+delta < 8000)) +
  geom_histogram(aes(x=bp), fill="darkorange", binwidth=10)  +
  geom_histogram(aes(x=ri), fill="darkorange", binwidth=10)  +
  theme(text = element_text(size=9, family="Helvetica"),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size=8, family="Helvetica"),
        axis.text.y = element_blank())
histo_L  

histo_M <- ggplot(data = del_df %>% filter(delta>1, seg=="M") %>% filter(length+delta < 8000)) +
  geom_histogram(aes(x=bp), fill="darkorange", binwidth=10)  +
  geom_histogram(aes(x=ri), fill="darkorange", binwidth=10)  +
  theme(text = element_text(size=9, family="Helvetica"),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size=8, family="Helvetica"),
        axis.text.y = element_blank())
histo_M  

histo_S <- ggplot(data = del_df %>% filter(delta>1, seg=="S") %>% filter(length+delta < 8000)) +
  geom_histogram(aes(x=bp), fill="darkorange", binwidth=10)  +
  geom_histogram(aes(x=ri), fill="darkorange", binwidth=10)  +
  theme(text = element_text(size=9, family="Helvetica"),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size=8, family="Helvetica"),
        axis.text.y = element_blank())
histo_S  

all_3_histo_p <- ggarrange(histo_L, histo_M, histo_S, ncol=3, nrow=1, widths=c(7,5,2))
all_3_histo_p
ggsave("dvg_histos.pdf", plot=all_3_histo_p, height=total_plot_height, width=(7*25.4), units="mm")
