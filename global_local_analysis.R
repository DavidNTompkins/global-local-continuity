##########
# Global-Local Preference Analysis
# David N. Tompkins

###### SECTION 0 - Libraries and Setup
# these lines will install the needed packages below
packages <- c("dplyr", "lme4", "emmeans", 
              "irrNA", "beepr", "ggplot2", "car", "lmerTest")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

library(dplyr)
library(tidyr)
library(lme4)
library(emmeans)
library(beepr)
library(ggplot2)
library(car)
library(lmerTest)
library(psych)
library(see)
library(ggrain)

###### SECTION 1 - READING IN DEIDENTIFIED DATA ######
# this is about 90mb of data - might take a  second to load.
main_path = gsub("global_local_analysis.R","",
                 rstudioapi::getSourceEditorContext()$path)

data_path = paste0(main_path, "/processed_data/")
aoi_path = paste0(main_path,"/aois/")

  rawdata<- read.csv(paste0(data_path,"gl_rawstream_data.csv")) # this reads in the dataset 
  aoi_set <- read.csv(paste0(aoi_path,"tight_fit.csv"),fileEncoding = "UTF-8-BOM") # reads in AOIs
  trial_index <- read.csv(paste0(aoi_path,"trial_index.csv"),fileEncoding = "UTF-8-BOM") # This provides an easy way to link looking between the first and second phases of a trial.
  

# Merging AOIs and then summarizing by trial. 
  rawdata2 <- left_join(rawdata,trial_index,by="stimulus") %>% left_join(.,aoi_set%>%select(-stimulus),by=c("trial_ordinal"))

  #summarizing by trial
  stim_summaries <- rawdata2 %>%
    filter(category_right !="Blink") %>% # Removing any blinks (uncategorized looking data '-' kept in, not limited to identified fixation and saccades)
    
    # Scoring each datapoint as in or out of each AOI
    rowwise()%>% # rowwise for aoi calcs - checking each sample for an aoi
    mutate(in_aoi1 = ifelse(between(xpos,aoi1_x, aoi1_x+aoi_w) & between(ypos,aoi1_y,aoi1_y+aoi_h),duration,0), # for clarity duration here is calcuated from timestamps but is very very close to 1/sample rate (1/120 ~= 8.33 ms)
           in_aoi2= ifelse(between(xpos,aoi2_x, aoi2_x+aoi_w) & between(ypos,aoi2_y,aoi2_y+aoi_h),duration,0)) %>% 
    ungroup() %>% # removing rowwise
    group_by(participant,stimulus) %>% # grouping to stim level
    reframe(dwell_left= sum(in_aoi1),
              dwell_right= sum(in_aoi2), # Summing total dwell time to each AOI
              global_side=global_side, # keeping these variables through summarizing - no change made
              local_side=local_side,
            trial_ordinal = trial_ordinal,
            trial_number=trial_number,
              stim_duration = Export.End.Trial.Time..ms., #length this stim was on screen
              trial_num = Trial, # kept primarily for ease of inspection
              tracking_ratio=Tracking.Ratio....) %>%
    distinct()
  



###### SECTION 3 - Applying Trial-level Exclusions ######

advance_record <- read.csv(paste0(data_path,"trial_advance_record.csv"),fileEncoding = "UTF-8-BOM") # This file contains a record of which trials were advanced by a look to the
# trigger AOI at runtime. A similar order can be created by examining looking to the 
# trigger AOI, but they are not precisely the same - possibly as the algo performed at runtime 
# does not perfectly match the algo used here (e.g., interpolation over samples). This file reflects the actual path participants took.
  
# Applying trial-level exclusions
  # Joining both records to data at stimulus level
  stim_summaries2 <-left_join(stim_summaries,advance_record, by =c("participant","stimulus")) #joining these files
  
  # cleaning read-in data
  stim_summaries2 <- stim_summaries2 %>%
    select(-X)%>%
    group_by(participant,trial_number) %>%
    mutate(advanced_by_fixation = ifelse(is.na(advanced_by_fixation),0,advanced_by_fixation), # swapping na values for secondary stim to 0 (i.e. we don't have information for these, the next line extends our info from sample stim for ease of use)
           advanced_by_fixation=max(advanced_by_fixation),
           ) # this  extends advanced_by_fixation from primary to secondary stim phases
  
  ## First filters out irrelevant lines (sample/first phase stim, and trials with no looking to either side), then uses runtime information to create exclusion
  stim_summaries5 <- stim_summaries2 %>%
    ungroup() %>%
    filter(grepl("Test",stimulus,fixed=TRUE)) %>% # removing non test-period data
    filter(advanced_by_fixation==1) %>%# Excluding trials without a contingent gaze
    filter(!(dwell_left==0 & dwell_right==0)) #cutting out trials with no look to either side.

  # creating exclusion stats:
  stim_summaries2 %>%     ungroup() %>%
    filter(grepl("Test",stimulus,fixed=TRUE)) %>% # removing non test-period data
    summarize(n_trials = n()) 
  
  stim_summaries2 %>%     ungroup() %>%
    filter(grepl("Test",stimulus,fixed=TRUE)) %>% # removing non test-period data
    filter(advanced_by_fixation==0) %>% #trials without a contingent gaze
  summarize(excluded_by_fixation = n())  
  
  stim_summaries2 %>%
    ungroup() %>%
    filter(grepl("Test",stimulus,fixed=TRUE)) %>% # removing non test-period data
    filter(advanced_by_fixation==1) %>%# Excluding trials without a contingent gaze
    filter((dwell_left==0 & dwell_right==0))%>% # trials with no looking
    summarize(no_left_right =n()) 
    
###### SECTION 4 - Creating Participant Summaries ######

# Reading in demographics
demographics <- read.csv(paste0(data_path,"gl_demographics.csv"),fileEncoding = "UTF-8-BOM")

# creating trial level summaries of left/right preference and global preference.
# As primary phase has been removed, here 1 stim is 1 trial
#  'ratio' is a misnomer for proportion (ratio of global to total looking, e.g., G/(L+G) )
  trial_summaries <- stim_summaries5 %>%
    group_by(participant,stimulus) %>%
    summarize(right_ratio=(dwell_right/(dwell_left+dwell_right)),
              global_ratio=ifelse(global_side=="RIGHT",right_ratio,(dwell_left/(dwell_left+dwell_right))),
              right_binary= round(right_ratio),
              global_binary=round(global_ratio),
              dwell_global = ifelse(global_side=="RIGHT",dwell_right,dwell_left),
              dwell_local = ifelse(global_side=="LEFT",dwell_right,dwell_left),
              dwell_left = dwell_left,
              dwell_right = dwell_right,
              trial_ordinal = trial_ordinal,
              global_side = global_side
    )%>%
    distinct()
  
  # Adding demographic information at trial level
  trial_summaries2 <- left_join(trial_summaries, demographics)
  
  # Creating participant summaries
  participant_summaries <- trial_summaries2 %>%
    group_by(participant)%>%
    reframe(global_preference= mean(global_ratio),
              right_preference = mean(right_ratio),
              ntrials = n(),
              sex = sex,
              age_months = age_months,
              age_group = age_group
    ) %>%
    distinct()
  
  mean(participant_summaries$ntrials)
  sd(participant_summaries$ntrials)
  range(participant_summaries$ntrials)
  
  number_of_participants<-nrow(participant_summaries)
  number_of_trials <- nrow(trial_summaries2)

# some sanity checks:
  participants <- trial_summaries2 %>%
    group_by(participant,age_group, age_months,sex,race,Hispanic,maternal_education) %>% 
    summarize(mean_gl = mean(global_ratio), #reasonable
              howmany_trials = n()) #1-20, yes
  
  # demos
  participants %>% ungroup() %>%
    summarise(
      n = n(),
      age_min = min(age_months, na.rm =T),
      age_max = max(age_months,na.rm=T),
      mean_age_months = mean(age_months, na.rm = TRUE),
      sd_age_months = sd(age_months, na.rm = TRUE),
      boys = sum(sex == "M", na.rm = TRUE),
      girls = sum(sex == "F", na.rm = TRUE)
    )
  
  participants %>%
    group_by(age_group) %>%
    summarise(
      n = n(),
      mean_age_months = mean(age_months, na.rm = TRUE),
      sd_age_months = sd(age_months, na.rm = TRUE),
      boys = sum(sex == "M", na.rm = TRUE),
      girls = sum(sex == "F", na.rm = TRUE)
    )
  
  participants %>% ungroup() %>%
    distinct(participant, race) %>%
    tidyr::separate_rows(race, sep = ",") %>%
    mutate(race = stringr::str_trim(race)) %>%
    distinct(participant, race) %>%
    count(race) %>%
    mutate(percent = round(n / n_distinct(participants$participant) * 100, 1))
  
  participants %>% ungroup() %>%
    select(Hispanic) %>%
    tidyr::pivot_longer(everything()) %>%
    count(name, value) %>%
    group_by(name) %>%
    mutate(percent = round(n / sum(n) * 100, 1))
  
  participants %>% ungroup() %>%
    select(maternal_education) %>%
    tidyr::pivot_longer(everything()) %>%
    count(name, value) %>%
    group_by(name) %>%
    mutate(percent = round(n / 80 * 100, 1)) # does not include the 2 unknown for base rates
  
  ###### SECTION 5 - Running analyses ######

  trial_summaries2$sex <- factor(trial_summaries2$sex,levels = c("M","F"))
  contrasts(trial_summaries2$sex) <- c(-.5,.5)
  

  # defining linear models
  m1 =  lmer(global_ratio ~ scale(age_months) + sex  + (1|participant)+ (1|stimulus), data=trial_summaries2)
  summary(m1)
  anova(m1)
  emmeans(m1, ~1) 
  emmeans(m1, ~1) |> test(null = 0.5)
  
  # second model removing participant
  m2 =  lmer(global_ratio ~ scale(age_months) + sex  + (1|stimulus), data=trial_summaries2)
  summary(m2)
  anova(m2)
  emmeans(m2, ~1) 
  emmeans(m2, ~1) |> test(null = 0.5)


  ##### Visualizing ###
  # Visualizing preference  (figure 2 in paper)

  age_group_labels <- c(12, 18, 24, 30, 36, 42, 48)
  overall_x <- length(age_group_labels) + 1 # for location of overall mean
  
  participant_summaries_plot <- participant_summaries %>%
    mutate(age_group_plot = as.numeric(factor(age_group, levels = age_group_labels)))
  
  overall_summary <- participant_summaries_plot %>%
    summarize(global_preference_se = sd(global_preference, na.rm = TRUE) / sqrt(sum(!is.na(global_preference))),
              global_preference = mean(global_preference, na.rm = TRUE),
              age_group_plot = overall_x)
  
  ggplot(data = participant_summaries_plot, aes(x = age_group_plot, y = global_preference)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40", linewidth = 0.7) +
    
    geom_boxplot(
      aes(x = age_group_plot - 0.14, group = age_group),
      width = 0.25, alpha = 0.3, outliers = FALSE,
      color = "darkred"
    ) +
    geom_jitter(aes(x = age_group_plot + 0.14, size=ntrials),color="#255F85", height = 0, width = 0.08, alpha = 0.5) +
    
    geom_pointrange(
      data = overall_summary,
      aes(x = age_group_plot, y = global_preference,
          ymin = global_preference - global_preference_se,
          ymax = global_preference + global_preference_se),
      inherit.aes = FALSE, size =0.6,
      color = "black", linewidth = 1.1
    ) +
    
    annotate("segment",
             x = c(0.5, overall_x - 0.3),
             xend = c(overall_x - 0.7, overall_x + 0.5),
             y = -Inf, yend = -Inf,
             linewidth = 0.5) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = c("0", ".25", ".50", ".75", "1.0"),
      name = "Global Preference (proportion)"
    ) +
    scale_x_continuous(
      limits = c(0.5, overall_x + 0.5),
      breaks = c(seq_along(age_group_labels), overall_x),
      labels = c(as.character(age_group_labels), "Overall"),
      name = "Age Group (months)"
    ) +
    scale_size(range = c(1, 4), name = "Participant\nTrials") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.line.x = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4)
    )
  
# Visualizing looking including looking to target.
rawdata3 <- rawdata2 %>%
  mutate(AOI_name = case_when(
    between(xpos,aoi1_x, aoi1_x+aoi_w) & between(ypos,aoi1_y,aoi1_y+aoi_h) ~ "LEFT",
    between(xpos,aoi2_x, aoi2_x+aoi_w) & between(ypos,aoi2_y,aoi2_y+aoi_h) ~ "RIGHT",
    between(xpos,target_x, target_x+aoi_w) & between(ypos,target_y,target_y+aoi_h) ~ "Target",
    
    .default = "Non-AOI"
  )) %>% mutate(
    cast_AOI = case_when(
      AOI_name == global_side ~ "Global",
      AOI_name == local_side ~ "Local",
      .default = AOI_name
    )
  )
  
# joining with stim_summaries 2 (last place we have familiarization and advancement trials)
advancement_record_details <- stim_summaries2 %>%ungroup() %>% select(participant, stimulus, advanced_by_fixation,trial_number)
rawdata4 <- rawdata3 %>% 
  left_join(., advancement_record_details, by=c("participant","trial_number","stimulus"))%>% filter(advanced_by_fixation==1) %>%
  ungroup() %>% group_by(participant,stimulus) %>% mutate(stim_start = min(timemark)) %>%ungroup() %>%
  group_by(participant,trial_number) %>% mutate(trial_start=max(stim_start)) %>% # basically, split trials by stim, grab the starts, then take the latter
  ungroup() %>% group_by(participant,trial_number) %>% mutate(time_from_trial_reveal=timemark-trial_start) %>%ungroup() %>%
  filter(time_from_trial_reveal<2001) %>% ungroup() %>%
  left_join(., demographics,by="participant") 

range(rawdata4$time_from_trial_reveal)

#simple time series approach
time_data <- rawdata4 %>%
  ungroup()  %>% mutate(bin_10 = floor(time_from_trial_reveal / 20) * 20 + 10) %>%
  group_by(bin_10) %>%
  summarize(
    looking_target = mean(cast_AOI == "Target"),
    looking_global = mean(cast_AOI == "Global"),
    looking_local = mean(cast_AOI == "Local")
  )


# figure 1 in paper
ggplot(data = time_data, aes(x = bin_10)) +
  geom_line(aes(y = looking_target, color = "Target"), linewidth = 1.3) +
  geom_line(aes(y = looking_global, color = "Global"), linewidth = 1.3) +
  geom_line(aes(y = looking_local, color = "Local"), linewidth = 1.3) +
  scale_color_manual(
    name ="AOI",
    values = c(
      "Target" = "blue",
      "Global" = "red",
      "Local" = "green"
    )
  ) +
  xlim(-500, 2000) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4)
  )+
  geom_vline(xintercept = 0)+
  labs(
    x = "Time relative to test phase start (ms)",
    y = "Proportion looking to AOI"
  )

## Does individual continued looking at target correlate with global or local preference?
# trim to test trials
# sum target looking prior to looking to global or local. Matching filtering done for original analysis
# within trials, does this relate to global preference

target_looking_trials <- rawdata4 %>%
  filter(grepl("Test",stimulus,fixed=TRUE)) %>% # test phases only
  group_by(participant,stimulus) %>% # grouping to trials
  filter(any(cast_AOI %in% c("Global","Local"))) %>% # only trials with local or global looking
  mutate(first_global_local=min(time_from_trial_reveal[cast_AOI %in% c("Global","Local")])) %>% #when did they look at the gl?
  reframe(target_looking=sum(ifelse(cast_AOI=="Target" & time_from_trial_reveal<first_global_local,duration,0)),
          ttgl = min(first_global_local)) %>%
  left_join(., trial_summaries2 %>%
              ungroup() %>%
              select(participant, stimulus, global_ratio,age_months,age_group,sex), # grabbing previously calculated global ratios
            by=c("participant","stimulus")) %>% filter(!is.na(sex)) # filtering out trials filtered in t_sums3


cor.test(target_looking_trials$ttgl,target_looking_trials$global_ratio) # not at all correlated - Note we report the correlation below, which is related to this (and inferrentially equiv.) as this violates iid

mean(target_looking_trials$ttgl) 
sd(target_looking_trials$ttgl) 
range(target_looking_trials$ttgl)

target_looking_participants <- target_looking_trials %>%
  ungroup() %>% group_by(participant,age_group,sex) %>%
  summarize(age_months=max(age_months),
            ttgl=mean(ttgl),
            mean_pref = mean(global_ratio))

cor.test(target_looking_participants$ttgl,target_looking_participants$age_months) # older kids faster to look 
cor.test(target_looking_participants$ttgl,target_looking_participants$mean_pref) # this is a more apt test of correlation between pref and looking to target (same inferential interp)

ggplot(data=target_looking_participants,aes(x=age_group,y=ttgl))+
  geom_boxplot(
    aes(group = age_group),
    width = 3, alpha = 0.3, outliers = FALSE,
    color = "gray30", fill = "gray80"
  ) +  geom_jitter(height=0,width=0.5,alpha=0.5,aes(color=sex))+
  geom_smooth(method="lm")+
  
  scale_color_manual(values = c("M" = "#7B2D8B", "F" = "#2E8B57"),
                     labels = c("M" = "Boys", "F" = "Girls"),
                     name = NULL) +
  scale_y_continuous(
    #limits = c(0, 1),
    #breaks = seq(0, 1, 0.25),
    #labels = c("0", ".25", ".50", ".75", "1.0"),
    name = "Latency to looking at global or local arrays (ms)"
  ) +
  scale_x_continuous(
    breaks = c(12, 18, 24, 30, 36, 42, 48),
    labels = c("12", "18", "24", "30", "36", "42", "48"),
    name = "Age Group (months)"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4)
  ) # not in paper, but kept for interest, note that visually similar to Figure 2, but different Y axis



#### Calibration Data Review #### 
# This concerns the calibration data collected during eye tracking. This also re-runs analyses after removing participants with errant or missing validation data.
# removing these does not affect results, which is clear below, but also implied by our singular fit with participant random effects.

cal_path <- paste0(main_path,"/calibration-data/")
validation <- read.csv(paste0(cal_path,"ValidationResults.txt"),sep = '\t')

# Validation includes originally entered participant IDs which have typos fixed in data above. Fixing here:
validation2 <- validation %>% 
  rename(participant = Participant) %>%
  mutate(participant = dplyr::recode(participant,
                              "P4083" = "DAVPIL4083",  # Fixing misentries of ids during testing 
                              "P4084" = "DAVPIL4084",
                              "DAVPIL4069a" = "DAVPIL4069",
                              "DAVPIL1043" = "DAVPIL4103")) %>%
    filter(participant !="DAVPIL4107_b", participant !="DAVPIL4017")%>% # excluding second file from a participant, as well as a participant with no data (e.g. did not begin task)
    mutate(participant = toupper(participant))%>%
    filter(participant != "DAVPIL4085", participant !="DAVPIL4036") #excluding based on color blindness

# aligning to included participants (i.e., removing 3 who had no contributed data)

validation3 <- full_join(participants,validation2,by="participant") %>% 
  mutate(right_x = as.numeric(Right.Eye.Deviation.X....),
         right_y = as.numeric(Right.Eye.Deviation.Y....),
         left_x = as.numeric(Left.Eye.Deviation.X....),
         left_y = as.numeric(Left.Eye.Deviation.Y....))

  checks <- validation3 %>% filter(is.na(Type)|is.na(sex)) # checking for matches. 2 have missing/ NA validation. One seemingly due to overwrite, the other is recorded as NA in 
  
mean(validation3$right_x, na.rm=T)
mean(validation3$right_y,na.rm=T)

sd(validation3$right_x, na.rm=T)
sd(validation3$right_y,na.rm=T)

validation4 <- validation3 %>%
  mutate(over_two =ifelse(right_x>=2.0 | right_y >=2.0 |is.na(right_x),1,0))

sum(validation4$over_two,na.rm=T) #2 na, 8 over 2. Unclear why these 8 were moved into testing

# recreating key analyses using only <2 degree validation error

exclude_validation <- validation4 %>% # only those we want to remove, making a list
  filter(over_two==1) %>%
  pull(participant)

trial_summaries_validation <- trial_summaries2 %>%
  filter(!participant %in% exclude_validation)

trial_summaries_validation$sex <- factor(trial_summaries_validation$sex,levels = c("M","F"))
contrasts(trial_summaries_validation$sex) <- c(-.5,.5)

n_distinct(trial_summaries_validation$participant)
nrow(trial_summaries_validation)

m1_validation = lmer(global_ratio ~ scale(age_months) + sex + (1|participant) + (1|stimulus),
                     data=trial_summaries_validation) # still singular
summary(m1_validation)
anova(m1_validation)
emmeans(m1_validation, ~1) |> test(null = 0.5)

m2_validation = lmer(global_ratio ~ scale(age_months) + sex + (1|stimulus),
                     data=trial_summaries_validation)
summary(m2_validation) # very much the same. 
anova(m2_validation)
emmeans(m2_validation, ~1)
emmeans(m2_validation, ~1) |> test(null = 0.5) # vs 0.55. both significant

target_looking_validation <- target_looking_participants %>%
  filter(!participant %in% exclude_validation)

cor.test(target_looking_validation$ttgl,target_looking_validation$age_months) # vs -0.51, both sig
cor.test(target_looking_validation$ttgl,target_looking_validation$mean_pref) # vs 0.002, both insig 

target_looking_trials_validation <- target_looking_trials %>%
  filter(!participant %in% exclude_validation)

mean(target_looking_trials_validation$ttgl) # vs 444.93
sd(target_looking_trials_validation$ttgl) # vs 259.53
         