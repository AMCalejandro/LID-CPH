

```{r}
###### GET THE VISIT NUMBER TO EXPLORE COGNITIVE DECLINE IN LID VS NON LID #######
###################################################################################
#saveRDS(timeDysk_updrsIII_Ldopadose %>% select(ID, visit_number), "../../GITHUB_REPOS/Dyskinesia_Project/data/Dyskinesias_visit_number.rds")  
lastvisit_VisitNumber_noDySkinesias = lastvisit_VisitNumber %>%
  filter(!ID %in% timeDysk_updrsIII_Ldopadose$ID) %>%
  rename(visit_number = Key) %>%
  mutate(visit_number = as.numeric(visit_number))
#saveRDS(timeDysk_updrsIII_Ldopadose %>% select(ID, visit_number), "../../GITHUB_REPOS/Dyskinesia_Project/data/Dyskinesias_visit_number.rds")  

PDLID_VISIT = bind_rows(lastvisit_VisitNumber_noDySkinesias, 
                        timeDysk_updrsIII_Ldopadose %>% 
                          select(ID, visit_number))
saveRDS(PDLID_VISIT, "../../GITHUB_REPOS/Dyskinesia_Project/data/PDLIDall_visitnumber.rds") 

###################################################################################
```

