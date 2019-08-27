library (stringr)
library (dplyr)
library (ggplot2)
library (gridExtra)
# library(cowplot)

# This script creates summaries of "split read" coverage and detected deletion DVGs for all genome segments
#
# Mark Stenglein, May 3, 2019

# import and tidy the data
if (!exists("depth_df")) source("./import_and_tidy_dvg_data.R")

# how may DVGs are there (by S, M, L)
# (just output to console)
del_df %>% group_by(seg) %>% summarize(n=n())

# this function creates labels for the y axis
# 
# It creates empty labels (" ") for any y axis value < 1
# This is because we are co-opting the Y-axis < 1 for plotting DVGs
yaxis_label_function <- function(breaks){
  labels <- vector(mode="character", length = 0)
  for (b in breaks){
    if (!is.na(b) & as.numeric(b) >= 1) {
      labels <- c(labels, b)
    }
    else
    {
      labels <- c(labels, " ")
    }
  }
  return (labels)
}

#
# main function to create plots showing coverage and DVGs
# 
# returns a ggplot 
combined_plot_function <- function(plot_depth_df, plot_del_df, logs_for_cartoons){
  
  # min/max DVGs plotted per genome
  max_dvg <- 60
  min_dvg <- 6
  
  # the number of logs of y axis to use to plot DVG cartoons
  if (missing(logs_for_cartoons)){
     logs_for_cartoons <- 2
  }
  max_log = 0
  min_log = max_log - logs_for_cartoons
  
  # these functions convert DVG# to a faux y axis log scale value
  # so that we can plot DVGs on the same y-axis as the coverage info
  row_to_logvalue <- function(row_num, n_dvg){
    log_break = abs(max_log - min_log) / n_dvg
    # dvg between 
    logval <- 10^(max_log-(log_break*row_num))
    return (logval)
  }
  
  logvalue_to_row <- function(logval, n_dvg){
    log_break = abs(max_log - min_log) / n_dvg
    row_num <- (max_log-log10(logval))/log_break
    return (row_num)
  }
  
  rect_delta <- function(logval, n_dvg){
    row_num <- logvalue_to_row(logval, n_dvg)
    return (row_to_logvalue(row_num + 0.25, n_dvg))
  }

  rect_delta_neg <- function(logval, n_dvg){
    row_num <- logvalue_to_row(logval, n_dvg)
    return (row_to_logvalue(row_num - 0.25, n_dvg))
  }
  
  # for debugging
  # plot_del_df <- filter(del_df, virus=="Fort_Sherman")
  
  # sort by supporting read counts and keep track of scaled value for plotting in a new column
  plot_del_df <- plot_del_df %>% 
    group_by (virus, seg) %>% 
    arrange(-counts, .by_group=TRUE) %>% 
    # mutate(n_dvg = max(min(max_dvg, n()), min_dvg), desc_order = row_to_logvalue(row_number(), n_dvg) )
    mutate(n_dvg = max(min(max_dvg, n()), min_dvg), row_num = row_number(), desc_order = row_to_logvalue(row_number(), n_dvg) ) %>%
    filter(row_num <= n_dvg)
  
  # main plot function
  p <- ggplot(plot_depth_df) +
    # total coverage
    geom_area(aes(x=position, y=depth), fill="grey90")  +
    # split reads coverage
    geom_area(aes(x=position, y=split_depth), fill="darkorange", alpha=0.9) + 
    # a white background for cartoons
    geom_rect(data=plot_del_df, aes(ymin=10^min_log, ymax=10^max_log, xmin=1, xmax=(length+delta)), fill="white") +
    # these geoms form the deleted genome cartoons 
    geom_segment(data=plot_del_df, 
                 aes(x=bp, xend=ri, y=desc_order, yend=desc_order, color=counts), 
                 alpha=1, size=0.3) +
    geom_rect(data=plot_del_df, 
              aes(xmin=0, xmax=bp, 
                  ymin=rect_delta_neg(desc_order, n_dvg), ymax=rect_delta(desc_order, n_dvg)), 
              color="black", fill="darkorange", size=0.1) +
    geom_rect(data=plot_del_df, 
              aes(xmin=ri, xmax=(length+delta), 
                  ymin=rect_delta_neg(desc_order, n_dvg), ymax=rect_delta(desc_order, n_dvg)), 
              color="black", fill="darkorange", size=0.1) +
    # this gradient colors the connecting lines for the deleted genomes
    scale_color_gradient(low='grey90',high='blue') +
    
    theme_bw() +
    theme(panel.grid = element_line(color="slategray1", linetype=5, size=0.15)) + 
    
    # this facets by virus (rows) and segment (co)
    facet_grid(virus~seg, scales="free", space="free_x") +
    scale_y_log10(limits=c(min_log,NA), labels=yaxis_label_function) +
    xlab("Genome position (nt)") +
    ylab("Coverage")
  return (p)
}

# This function will create plots for the specified viruses
plot_viruses <- function(virus_names, logs_for_cartoons){
  
  # subset the main dataframes to get the data just for these viruses
  subset_di <- del_df %>% filter(virus %in% virus_names)
  subset_cov <- combined_depth_df %>% filter(virus %in% virus_names)
  
  # output the virus names to the console
  print(paste0(virus_names))
  
  # p returned here is a ggplot object
  p <- combined_plot_function(subset_cov, subset_di, logs_for_cartoons)
  
  # print p will make the plots appear on viewer
  print(p)
}

# test some plotting: note this can be slow
plot_viruses(c("Bunyamwera", "Main_Drain_72"))
plot_viruses(c("Yaba_7", "Turlock", "Moriche", "Turlock"))
plot_viruses(c("Yaba_7", "San_Juan"))
plot_viruses(c("Yaba_7", "Moriche"))
plot_viruses(c("Yaba_7", "Fort_Sherman"))
plot_viruses(c("Anadyr", "Birao", "Bujaru", "San_Juan", "Yaba_7", "Moriche"))
plot_viruses(c("Yaba_7"))
plot_viruses(c("Moriche"))
plot_viruses(c("Cananeia", "Mirim"))

# This function will break up a bunch of viruses into smaller subsets for plot generation
# the reason for doing this is to create multi-page PDF output
plot_some_viruses <- function(viruses, plots_per_page, pdf_name, logs_for_cartoons){

  # open a PDF for writing
  if (!missing(pdf_name)) {
    pdf(pdf_name, 7, 10)
  }
  
  if (missing(logs_for_cartoons)){
    logs_for_cartoons = 2
  }
  
  # iterate through the viruses, doing so many per page
  for (i in seq(1, length(viruses), plots_per_page)) {
    
    # plots_per_page at a time
    subset_viruses <- viruses[i:(i+(plots_per_page-1))]
    
    # output to console which ones we're doing
    print(paste0(subset_viruses))
    
    # generate & print plot 
    plot_viruses(subset_viruses, logs_for_cartoons)
  }
  if (!missing(pdf_name)) {
    dev.off()
  }
  # 
}

# test the function --> should produce a pdf named can_mir.pdf
plot_some_viruses(c("Cananeia", "Mirim"), 8, "can_mir.pdf", 2)

#
# This will create plots for all viruses
#
# list of segments w/ DI-tect DVGs
viruses <- combined_depth_df %>% group_by(virus) %>% summarize()
plot_some_viruses(viruses$virus, 9, "all_viruses_dvg_plot.pdf", 3)

# Manually list viruses with evidence of DVGs for making figures
#
# In some cases, add n extra segments just so individual plots/grids have similar heights
# will remove them manually in Adobe Illustrator later.
#
# Plot those viruses with more DVGs and fewer DVGs separtely - give those with 
# more more y axis space for the DVG cartoons.
#
# Also plot nairoviruses separately from orthobunyaviruses since their L segments are so long 

# plot_some_viruses <- function(viruses, plots_per_page, pdf_name, logs_for_cartoons){
# Those viruses that exhibited many DVG forms
main_drain_72 <- c("Main_Drain_72", "Bozo", "Fort_Sherman", "Kaisodi", "Moriche", "M_Poko") 
plot_some_viruses(main_drain_72, 6, "main_drain_dvg_plot.pdf", 4)

more_dvg_forms <- c("Anadyr", "Bozo", "Fort_Sherman", "Kaisodi", "Moriche", "M_Poko", "Oriboca", "San_Juan", "Showke", "Tensaw", "Vinces", "Yaba_7", "Birao")
# more_dvg_forms <- c("Anadyr", "Bozo", "Fort_Sherman", "Kaisodi", "Main_Drain_72", "Moriche", "M_Poko", "Oriboca", "San_Juan", "Showke", "Tensaw", "Vinces")
plot_some_viruses(more_dvg_forms, 6, "more_dvg_plot.pdf", 4)

non_nairo_dvgs <- c("Anadyr", "Anhanga", "Anopheles_B", "Bakau", "Belem", "Benefica", 
                "Birao", "Boraceia_3", "Bozo", "Bujaru",  "Bunyamwera", "Caimito", 
                "Cananeia", "Caraparu", "Fort_Sherman", "Kaisodi", "Lednice", 
                "Lokern", "Main_Drain", "Yaba_7", "Marituba", "Mirim",  "Moriche", "M_Poko", "Oriboca", 
                "Palma", "Potosi", "Punta_Toro",  "San_Juan", "Showke", 
                "Tacaiuma", "Tensaw", "Timboteua", "Tlacotalpan", "Turlock", "Vinces", 
                "Weldona", "Wongal", "Wyeomia", "Yaba_7")
non_nairo_dvgs_md <- c("Anadyr", "Anhanga", "Anopheles_B", "Bakau", "Belem", "Benefica", 
                "Birao", "Boraceia_3", "Bozo", "Bujaru",  "Bunyamwera", "Caimito", 
                "Cananeia", "Caraparu", "Fort_Sherman", "Kaisodi", "Lednice", 
                "Lokern", "Main_Drain", "Main_Drain_72", "Marituba", "Mirim",  "Moriche", "M_Poko", "Oriboca", 
                "Palma", "Potosi", "Punta_Toro",  "San_Juan", "Showke", 
                "Tacaiuma", "Tensaw", "Timboteua", "Tlacotalpan", "Turlock", "Vinces", 
                "Weldona", "Wongal", "Wyeomia", "Yaba_7")
plot_some_viruses(non_nairo_dvgs, 8, "non_nairo_dvg_plot.pdf", 3)
plot_some_viruses(non_nairo_dvgs_md, 8, "non_nairo_dvg_plot_md.pdf", 3)

# Nairos - long L segments
# put in extra segments just so individual plots/grids have similar heights
# will remove them manually in AI later
nairo_dvgs <- c("Hughes", "Yogue", "Anhanga", "Anopheles_B", "Bakau", "Belem")
plot_some_viruses(nairo_dvgs, 6, "nairo_dvg_plot.pdf", 1)

# fewer DVG forms
fewer_dvgs <- c("Anhanga", "Anopheles_B", "Bakau", "Belem", "Benefica", "Birao", "Boraceia_3", "Bujaru",  "Bunyamwera", "Caimito", 
                "Cananeia", "Caraparu", "Lednice", "Lokern", "Main_Drain", "Marituba", "Mirim",  "Palma", "Potosi", "Punta_Toro", 
                "Tacaiuma", "Timboteua", "Tlacotalpan", "Turlock", "Weldona", "Wongal", "Wyeomia", "Yaba_7")
plot_some_viruses(fewer_dvgs, 10, "fewer_dvg_plot.pdf", 2)

