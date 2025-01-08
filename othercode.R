jmex = wiDB_data(projects = "00044", types = "Tap")
jmexd = jmex$data

write.csv(jmexd, "Mexicotapcj.csv")

summproject <- twd.keep %>%
  group_by(Project_ID)%>%
  summarize(n=n())


print(tw$project_ids)
list(tw)
list(tw.keep)
