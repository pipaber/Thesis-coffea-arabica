library(gt)
library(tidyverse)
library(janitor)
library(viridis)

table_filt_reads_1<-fastp_filtered_reads_plot%>%
  clean_names()%>%
    summarize_at(.vars = c("passed_filter",
                         "low_quality",
                         "too_many_n",
                         "too_short"),
               .funs = list(Promedio=mean))%>%
  mutate(funcion="Promedio")
  
table_filt_reads_2<-fastp_filtered_reads_plot%>%
  clean_names()%>%
  summarize_at(.vars = c("passed_filter",
                         "low_quality",
                         "too_many_n",
                         "too_short"),
               .funs = list(Promedio=median))%>%
  mutate(funcion="Mediana")
  
table_join_filt<-
  table_filt_reads_1%>%
  bind_rows(table_filt_reads_2)
  
table_join_filt%>%
  gt(rowname_col = "funcion")%>%
  tab_header(title = md("**<u>Tabla 1:</u> Estadísticas resumen de *reads***"))%>%
  opt_align_table_header(align = "left")%>%
  fmt_scientific(
    columns = everything()  
  )%>%
  opt_table_font(font = "TimesNewRoman")%>%
  cols_label(
    passed_filter_Promedio = md("Pasaron filtro"),
    low_quality_Promedio = md("Baja calidad"),
    too_many_n_Promedio = md("Muchos N"),
    too_short_Promedio = md("Muy cortos")
    )%>%
  gtsave("tab1.rtf")  

fastp_seq_quality_plot%>%
  pivot_longer(SRR11711678:SRR11711951,names_to = "Muestras",values_to = "Calidad")%>%
  ggplot(aes(`Read Position`,Calidad,group=Muestras,colour=Muestras))+
  geom_line()+
  coord_cartesian(ylim=c(0,50))+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(x="Posición de las bases",y="Calidad (%)")

fastp_filtered_reads_plot %>%
  pivot_longer("Passed Filter":"Too short",
               names_to = "Categoria",
               values_to = "Value") %>%
  arrange(Categoria, desc(Value)) %>%
  ggplot(aes(y = Category, fill = Categoria, x = Value)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(
    discrete = T,
    labels = c("Baja calidad", "Pasaron el filtro", "Muchos N", "Muy cortos")
  )+
  labs(fill="Estadísticas de filtrado",y="Muestras",x="Número de *reads*")+
  theme(
    axis.title.x = ggtext::element_markdown(),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "gray90"),
    panel.grid.minor = element_line(colour = "gray90")
  )
  


table_qc_1<-
  general_stats_qc%>%
  clean_names()%>%
  summarize_at(.vars=c("percent_duplication","percent_q30",
                       "gc_content", "percent_pf"),
               .funs= list(Promedio=mean))%>%
  mutate(funcion = "Promedio")


table_qc_2<-
  general_stats_qc%>%
  clean_names()%>%
  summarize_at(.vars=c("percent_duplication","percent_q30",
                       "gc_content", "percent_pf"),
               .funs= list(Promedio=median))%>%
  mutate(funcion = "Mediana")

table_join_qc<-
  table_qc_1%>%
  bind_rows(table_qc_2)

table_join_qc%>%
  gt(rowname_col = "funcion")%>%
  tab_header(title = md("**<u>Tabla 2</u>: Estadísticas resumen de la calidad de *reads***"))%>%
  opt_align_table_header(align = "left")%>%
  fmt_percent(columns = everything())%>%
  opt_table_font(font = "TimesNewRoman")%>%
  cols_label(
    percent_duplication_Promedio = md("Porcentaje de duplicación"),
    percent_q30_Promedio = md("Porcentaje > Q30"),
    gc_content_Promedio = md("Porcentaje de GC"),
    percent_pf_Promedio = md("Porcentaje de reads que<br>pasaron el filtro")
  )%>%gtsave("tabqc.rtf")
  

