library(here)
library(data.table)
library(ggplot2)

comp <- fread(here("rawdata/ADCONSP2_2023-01-15_06-27-54.csv"))
setkey(comp, `Participant ID`)

comp <- comp[`Study completion status (char)` == "Completed"]

compsum <- comp[, .(pid = `Participant ID`,
                    trt = `Treatment Group (Char)`, 
                    stratum = Stratum,
                    rand_date = `Date of Randomization`, 
                    disc_date = `Date of Study Intervention Discont.`,
                    lastfu_date = `Date of Last Follow-Up`, 
                    age = `Age (months)`,
                    sex = `Sex (char)`,
                    race = `Primary Ethnicity`,
                    itt = `Intent-to-treat (ITT) Sample`, 
                    fu_days = `Follow up Time (days)`,
                    grams_per_week = `Grams Consumed Per Week`
                    )]

ggplot(compsum, aes(x = grams_per_week, fill = trt)) + 
  geom_histogram(position = "dodge") 

ggplot(compsum, aes(x = grams_per_week, fill = trt)) + geom_histogram(position = "dodge")+ theme_bw() + 
  theme(legend.position = "bottom") + 
  scale_fill_discrete("Assigned treatment") +
  scale_x_continuous("Grams of peanuts consumed per week")
ggsave(here("results/peanut-hist.pdf"), width = 6.5, height = 4.25)


by(compsum$grams_per_week, compsum$trt, summary)

compsum[, consume_group := cut(compsum$grams_per_week, c(-1, .2, 6, 10, Inf), 
                               labels= c("less than .2g", ".2 to 6", "6 to 10", "10 or more"))]


table(compsum$consume_group, compsum$trt)

## outcome

outc <- fread(here("rawdata/ADPOUT1_2020-06-25_07-59-05.csv"))
#opp <- fread(here("rawdata/ADPP2_2020-06-25_13-01-37.csv"))

data <- merge(compsum, outc[, .(pid = `Participant ID`, 
         reaction = `Reaction`, 
         outcome= `Outcome of Peanut OFC at V60 (character)` == "Allergic", 
         outcomeraw= `Outcome of Peanut OFC at V60 (character)`)], 
      by = "pid", all.x = TRUE, all.y = FALSE)[outcomeraw %in% c("Allergic", "Tolerant")]

## old definition of treatment taken


rawpeanut <- read.csv2(here("rawdata/ADPP2_2020-06-25_13-01-37.csv"))

rawpeanut$X2 <- ifelse(rawpeanut$Treatment.Group..Char. == "Peanut Consumption",
             rawpeanut$Num..times.consumed.more.than.3.grams. != "0.0",
             rawpeanut$Numb..of.times..gt..0.2.grams..avoiders. != "0.0")

data <- merge(data, rawpeanut[, c("Participant.ID", "X2")],
      by.x="pid", by.y = "Participant.ID", all.x = TRUE, all.y = FALSE)

saveRDS(data, file = here("data/analysis_data.rds"))
