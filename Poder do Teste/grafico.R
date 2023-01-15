
library(tidyverse)
library(patchwork)
library(ggpubr)

# Sem Bootstrap ####

load('C:/Users/Pedro/Desktop/Pedro/UFG/Monografia/Codigos/power/trv_power_sb.RData')
load('C:/Users/Pedro/Desktop/Pedro/UFG/Monografia/Codigos/power/wald_power_sb.RData')
load('C:/Users/Pedro/Desktop/Pedro/UFG/Monografia/Codigos/power/score_power_sb.RData')
load('C:/Users/Pedro/Desktop/Pedro/UFG/Monografia/Codigos/power/grad_power_sb.RData')

# Com Bootstrap ####

load('C:/Users/Pedro/Desktop/Pedro/UFG/Monografia/Codigos/power/trv_power_cb.RData')
load('C:/Users/Pedro/Desktop/Pedro/UFG/Monografia/Codigos/power/wald_power_cb.RData')
load('C:/Users/Pedro/Desktop/Pedro/UFG/Monografia/Codigos/power/score_power_cb.RData')
load('C:/Users/Pedro/Desktop/Pedro/UFG/Monografia/Codigos/power/grad_power_cb.RData')


# Teste Razão Verossimilhança ####


trv_power_cb <- bind_cols(
  trv_n10
  , trv_n20 %>% .$value
  , trv_n30 %>% .$value
  , trv_n50 %>% .$value
) %>% 
  magrittr::set_colnames(
    c(
      'parametro'
      ,'n = 10'
      ,'n = 20'
      ,'n = 30'
      ,'n = 50'
    )
  ) %>% 
  reshape2::melt(id.vars = 'parametro') %>% 
  ggplot(aes(x = parametro, y = value, color = variable))+
  geom_line()+geom_point()+
  labs(
    title="Com Bootstrap"
    , x = "Probabilidade de Rejeição"
    , y = " "
    , color = 'Tamanho amostral'
  )+
  scale_colour_grey()+
  theme_light()+
  theme(
    legend.position = 'bottom'
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , plot.title = element_text(hjust = 0.5)
  )

trv_power_sb <- trv_power_sb+
  scale_colour_grey()+
  labs(
    title="Sem Bootstrap"
    , x = "Probabilidade de Rejeição"
    , y = "Valor do Parâmetro"
    , color = 'Tamanho amostral'
    , legend.grob = 'NULL'
  )

p1 <- ggarrange(
  trv_power_sb
  , trv_power_cb
  , ncol = 2
  , nrow = 1
  , common.legend = TRUE
  , legend = "bottom"
)

p1 <- annotate_figure(
  p1
  , top = text_grob(
    "Teste da razão de verossimilhança"
    , color = "black"
      , face = "bold"
    , size = 14
  )
)


# Teste de Wald ####


wald_power_cb <- bind_cols(
  wald_n10
  , wald_n20 %>% .$value
  , wald_n30 %>% .$value
  , wald_n50 %>% .$value
) %>% 
  magrittr::set_colnames(
    c(
      'parametro'
      ,'n = 10'
      ,'n = 20'
      ,'n = 30'
      ,'n = 50'
    )
  ) %>% 
  reshape2::melt(id.vars = 'parametro') %>% 
  ggplot(aes(x = parametro, y = value, color = variable))+
  geom_line()+geom_point()+
  labs(
    title="Com Bootstrap"
    , x = "Probabilidade de Rejeição"
    , y = " "
    , color = 'Tamanho amostral'
  )+
  scale_colour_grey()+
  theme_light()+
  theme(
    legend.position = 'bottom'
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , plot.title = element_text(hjust = 0.5)
  )

wald_power_sb <- wald_power_sb+
  scale_colour_grey()+
  labs(
    title="Sem Bootstrap"
    , x = "Probabilidade de Rejeição"
    , y = "Valor do Parâmetro"
    , color = 'Tamanho amostral'
    , legend.grob = 'NULL'
  )

p2 <- ggarrange(
  wald_power_sb
  , wald_power_cb
  , ncol = 2
  , nrow = 1
  , common.legend = TRUE
  , legend = "bottom"
)

p2 <- annotate_figure(
  p2
  , top = text_grob(
    "Teste de Wald"
    , color = "black"
      , face = "bold"
    , size = 14
  )
)


# Teste de Score ####


score_power_cb <- bind_cols(
  score_n10
  , score_n20 %>% .$value
  , score_n30 %>% .$value
  , score_n50 %>% .$value
) %>% 
  magrittr::set_colnames(
    c(
      'parametro'
      ,'n = 10'
      ,'n = 20'
      ,'n = 30'
      ,'n = 50'
    )
  ) %>% 
  reshape2::melt(id.vars = 'parametro') %>% 
  ggplot(aes(x = parametro, y = value, color = variable))+
  geom_line()+geom_point()+
  labs(
    title="Sem Bootstrap"
    , x = "Probabilidade de Rejeição"
    , y = "Valor do Parâmetro"
    , color = 'Tamanho amostral'
  )+
  scale_colour_grey()+
  theme_light()+
  theme(
    legend.position = 'bottom'
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , plot.title = element_text(hjust = 0.5)
  )

score_power_sb <- score_power_sb+
  scale_colour_grey()+
  labs(
    title="Com Bootstrap"
    , x = "Probabilidade de Rejeição"
    , y = "Valor do Parâmetro"
    , color = 'Tamanho amostral'
    , legend.grob = 'NULL'
  )

p3 <- ggarrange(
  score_power_sb
  , score_power_cb
  , ncol = 2
  , nrow = 1
  , common.legend = TRUE
  , legend = "bottom"
)

p3 <- annotate_figure(
  p3
  , top = text_grob(
    "Teste Escore"
    , color = "black"
      , face = "bold"
    , size = 14
  )
)


# Teste Gradiente ####


grad_power_cb <- bind_cols(
  grad_n10
  , grad_n20 %>% .$value
  , grad_n30 %>% .$value
  , grad_n50 %>% .$value
) %>% 
  magrittr::set_colnames(
    c(
      'parametro'
      ,'n = 10'
      ,'n = 20'
      ,'n = 30'
      ,'n = 50'
    )
  ) %>% 
  reshape2::melt(id.vars = 'parametro') %>% 
  ggplot(aes(x = parametro, y = value, color = variable))+
  geom_line()+geom_point()+
  labs(
    title="Com Bootstrap"
    , x = "Probabilidade de Rejeição"
    , y = "Valor do Parâmetro"
    , color = 'Tamanho amostral'
  )+
  scale_colour_grey()+
  theme_light()+
  theme(
    legend.position = 'bottom'
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , plot.title = element_text(hjust = 0.5)
  )

grad_power_sb <- grad_power_sb+
  scale_colour_grey()+
  labs(
    title="Sem Bootstrap"
    , x = "Probabilidade de Rejeição"
    , y = "Valor do Parâmetro"
    , color = 'Tamanho amostral'
    , legend.grob = 'NULL'
  )

p4 <- ggarrange(
  grad_power_sb
  , grad_power_cb
  , ncol = 2
  , nrow = 1
  , common.legend = TRUE
  , legend = "bottom"
)

p4 <- annotate_figure(
  p4
  , top = text_grob(
    "Teste Gradiente"
    , color = "black"
      , face = "bold"
    , size = 14
  )
)



p1
p2
p3
p4


