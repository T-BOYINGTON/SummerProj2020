
```{r}
# library
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(splitstackshape)
library(vegan)
library(ggeffects)
# library()
```

```{r}
setwd("~/Documents/SummerProject2020")
```


///////////////////////////Below == sel coef 0.014///////////////////////////


```{r}
mut2<- read.table ('build/test_bottleneck_p1_mutSummary_m2.txt', header = TRUE)
# ^^ deletrious mutation
mut3<- read.table ('build/test_bottleneck_p1_mutSummary_m3.txt', header = TRUE)
# ^^beneficial mutation
mut1<- read.table ('build/test_bottleneck_p1_mutSummary_m1.txt', header = TRUE) 
# ^^Natural mutation

p1ind<- read.table ('build/test_bottleneck_p1_ind_fitness.txt', header = TRUE)
p1het <- p1ind
```


```{r}
gen.pop <- p1ind %>%
select(gen, popSize)
 gen.pop1 <- ggplot(gen.pop, aes(x = gen, y = popSize)) +
          geom_jitter(width = 0.1, height = 0) + labs()
gen.pop1

# population size over generations
```

```{r}
dfm2 <- mut2 %>%
  select(gen, freqMut)

dfm2 <- dfm2 %>%   cSplit(dfm2, splitCols = "freqMut", sep = '_', direction = "long", fixed = TRUE, drop = FALSE, makeEqual= TRUE, type.convert = TRUE)
# splits freqMut col wherever a '_' occurs, allowing study of frequencies to be made
p1ind <- p1ind %>%   cSplit(p1ind, splitCols = "indvFitness", sep = '_', direction = "long", fixed = TRUE, drop = FALSE, makeEqual= TRUE, type.convert = TRUE)

p1het <- p1het %>%
  cSplit(p1het, splitCols = "Indv_Ho", sep = '_', direction = "long", fixed = TRUE, drop = FALSE, makeEqual= TRUE, type.convert = TRUE) 
p1het<-p1het %>%   mutate(het = (1-p1het$Indv_Ho))
p1het<-p1het %>% select(gen, Indv_Ho, het) %>% 
rename("Hom" = "Indv_Ho")
p1ind <- bind_cols(p1ind, p1het) %>% 
   select(!("gen...11")) %>% 
    rename("gen" = "gen...1")
remove(p1het)

dfm2<- dfm2 %>% 
  mutate(d.category = case_when(
    freqMut >= 0 & freqMut <= 0.3 ~ "Low",
    freqMut >0.3 & freqMut <= 0.6 ~ "Medium",
    freqMut >0.6 & freqMut < 1.0 ~ "High",
    freqMut == 1.0 ~ "Fixed"
  ))
# categorises m2 frequency into level of effect they may have

gen.freqM1 <- ggplot(mut1, aes(x = gen, y = Mutcount)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Mutation count (M1)") + theme_classic()

gen.freqM2 <- ggplot(mut2, aes(x = gen, y = Mutcount)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Deletrious Mutation count (M2)") +theme_classic()

gen.freqM3 <- ggplot(mut3, aes(x = gen, y = Mutcount)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Beneficial Mutation count (M3)") +theme_classic()

gen.freqM1
gen.freqM2
gen.freqM3
  # plots freq of each mut type. m1 and m3 show common upward trend while m2 is more rollercoastered (any selection would remove del mut whereas m1 and m3 are either neutral or actively beneficial) 
```


```{r}
dfm2 <- dfm2 %>% 
  mutate(freqMut = as.numeric(freqMut)) %>% 
filter(!is.na(freqMut)) %>% 
       mutate  (indvL = -((-0.014)*(2*0.5* freqMut +(1-(2*0.5)*((freqMut)^2))))
)   
# ^^ works out individual genetic load; 
# below works out additive load vv
# -(sum(del_selCoef*(2*del_domCoef*fMut+(1-(2*del_domCoef)*((fMut)^2)))))
dfm2 <- dfm2 %>% 
group_by(gen) %>% 
mutate(sum_L_by_gen = (sum(indvL)))

m2fixed<- dfm2 %>%
  filter(freqMut == 1.0)
#  fixed where freqMut == 1.0,
m2seg<- dfm2 %>% 
  filter(!(freqMut == 1.0))
#  otherwise segregating.

Sum.D.L <- ggplot(dfm2, aes(x = gen, y = sum_L_by_gen, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs() +theme_classic()
D.L <- ggplot(dfm2, aes(x = gen, y = indvL, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs() +theme_classic()
D.L
Sum.D.L  
# plot of indv genetic load /generation (D.L) and additive load/generation (Sum.D.L)

m2Fload.gen <- (ggplot(m2fixed, aes(x= gen, y= sum_L_by_gen, colour= freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Sum Fixed Load") +theme_classic())

m2Sload.gen <- (ggplot(m2seg, aes(x = gen, y = sum_L_by_gen, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Sum Segregation Load") +theme_classic())

m2Fload.gen  # if blank then no fixed M2 mutations
# Drift Load/generation
m2Sload.gen
# Segregation Load/generation

# all colour coded by allele frequency
```

```{r}

tdfm2 <- dfm2 %>% 
group_by("d.category") %>% 
count(d.category)
# counts n of each category ("Low", "Medium", and "High")
p1m2.freq <- (ggplot(dfm2, aes(x= gen, y= freqMut, colour= d.category)) +
          geom_jitter  (width = 0.1, height = 0) + labs( y = "Mutation frequency") +theme_classic())
p1m2.freq

p1count <- dfm2 %>% 
group_by(gen) %>% 
count(d.category)

p1m2.count<- ggplot(p1count, aes(x = gen, y = n, colour = d.category)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Deletrious Mutation count (M2)") +theme_classic()
p1m2.count

# plots the n of categories per gen

```

/////////////Below == sel coef 0.05/////////////

```{r}
p2mut2<- read.table ('build/test_bottleneck_p2_mutSummary_m2.txt', header = TRUE)
# ^^ deletrious mutation
p2mut3<- read.table ('build/test_bottleneck_p2_mutSummary_m3.txt', header = TRUE)
# ^^beneficial mutation
p2mut1<- read.table ('build/test_bottleneck_p2_mutSummary_m1.txt', header = TRUE) 
# ^^Natural mutation

p2ind<- read.table ('build/test_bottleneck_p2_ind_fitness.txt', header = TRUE)

p2het <- p2ind

gen.pop2 <- p2ind %>%
select(gen, popSize)
 gen.pop2 <- ggplot(gen.pop2, aes(x = gen, y = popSize)) +
          geom_jitter(width = 0.1, height = 0) + labs()
gen.pop2

```

```{r}
p2dfm2 <- p2mut2 %>%
  select(gen, freqMut)

p2dfm2 <- p2dfm2 %>%   cSplit(p2dfm2, splitCols = "freqMut", sep = '_', direction = "long", fixed = TRUE, drop = FALSE, makeEqual= TRUE, type.convert = TRUE)
# splits freqMut col wherever a '_' occurs, allowing study of frequencies to be made
p2ind <- p2ind %>%   cSplit(p2ind, splitCols = "indvFitness", sep = '_', direction = "long", fixed = TRUE, drop = FALSE, makeEqual= TRUE, type.convert = TRUE)

p2het <- p2het %>%
  cSplit(p2het, splitCols = "Indv_Ho", sep = '_', direction = "long", fixed = TRUE, drop = FALSE, makeEqual= TRUE, type.convert = TRUE) 
p2het<-p2het %>%   mutate(het = (1-p2het$Indv_Ho))
p2het<-p2het %>% select(gen, Indv_Ho, het) %>% 
rename("Hom" = "Indv_Ho")
p2ind <- bind_cols(p2ind, p2het) %>% 
   select(!("gen...11")) %>% 
    rename("gen" = "gen...1")
remove(p2het)

p2dfm2<- p2dfm2 %>% 
  mutate(d.category = case_when(
    freqMut >= 0 & freqMut <= 0.3 ~ "Low",
    freqMut >0.3 & freqMut <= 0.6 ~ "Medium",
    freqMut >0.6 & freqMut < 1.0 ~ "High",
    freqMut == 1.0 ~ "Fixed"
  ))
# categorises m2 frequency into level of effect they may have

p2gen.freqM1 <- ggplot(p2mut1, aes(x = gen, y = Mutcount)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Mutation count (M1)") + theme_classic()

p2gen.freqM2 <- ggplot(p2mut2, aes(x = gen, y = Mutcount)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Deletrious Mutation count (M2)") +theme_classic()

p2gen.freqM3 <- ggplot(p2mut3, aes(x = gen, y = Mutcount)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Beneficial Mutation count (M3)") +theme_classic()

p2gen.freqM1
p2gen.freqM2
p2gen.freqM3
  # plots freq of each mut type. m1 and m3 show common upward trend while m2 is more rollercoastered (any selection would remove del mut whereas m1 and m3 are either neutral or actively beneficial) 
```


```{r}
p2dfm2 <- p2dfm2 %>% 
  mutate(freqMut = as.numeric(freqMut)) %>% 
filter(!is.na(freqMut)) %>% 
       mutate  (indvL = -((-0.05)*(2*0.5* freqMut +(1-(2*0.5)*((freqMut)^2))))
)   
# ^^ works out individual genetic load; 
# below works out additive load vv
# -(sum(del_selCoef*(2*del_domCoef*fMut+(1-(2*del_domCoef)*((fMut)^2)))))
p2dfm2 <- p2dfm2 %>% 
group_by(gen) %>% 
mutate(sum_L_by_gen = (sum(indvL)))

p2m2fixed<- p2dfm2 %>%
  filter(freqMut == 1.0)
#  fixed where freqMut == 1.0,
p2m2seg<- p2dfm2 %>% 
  filter(!(freqMut == 1.0))
#  otherwise segregating.

p2Sum.D.L <- ggplot(p2dfm2, aes(x = gen, y = sum_L_by_gen, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs() +theme_classic()
p2D.L <- ggplot(p2dfm2, aes(x = gen, y = indvL, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs() +theme_classic()
p2D.L
p2Sum.D.L  # if blank then no fixed M2 mutations
# plot of indv genetic load /generation (D.L) and additive load/generation (Sum.D.L)

p2m2Fload.gen <- (ggplot(p2m2fixed, aes(x = gen, y = sum_L_by_gen, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Sum Fixed Load") +theme_classic())

p2m2Sload.gen <- (ggplot(p2m2seg, aes(x = gen, y = sum_L_by_gen, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Sum Segregation Load") +theme_classic())

p2m2Fload.gen 
# Drift Load/generation
p2m2Sload.gen
# Segregation Load/generation

```

```{r}
p2tdfm2 <- p2dfm2 %>% 
group_by("d.category") %>% 
count(d.category)
# counts n of each category ("Low", "Medium", and "High")
p2m2.freq <- (ggplot(p2dfm2, aes(x= gen, y= freqMut, colour= d.category)) +
          geom_jitter  (width = 0.1, height = 0) + labs( y = "mutation frequency") +theme_classic())
p2m2.freq

p2count <- p2dfm2 %>% 
group_by(gen) %>% 
count(d.category)

p2m2.count<- ggplot(p2count, aes(x = gen, y = n, colour = d.category)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Deletrious Mutation count (M2)") +theme_classic()
p2m2.count
```


///////////////Below == sel coef is 0.01////////////////////

```{r}
p3mut2<- read.table ('build/test_bottleneck_p3_mutSummary_m2.txt', header = TRUE)
# ^^ deletrious mutation
p3mut3<- read.table ('build/test_bottleneck_p3_mutSummary_m3.txt', header = TRUE)
# ^^beneficial mutation
p3mut1<- read.table ('build/test_bottleneck_p3_mutSummary_m1.txt', header = TRUE) 
# ^^Natural mutation

p3ind<- read.table ('build/test_bottleneck_p3_ind_fitness.txt', header = TRUE)
p3het <- p3ind
```


```{r}
gen.pop3 <- p3ind %>%
select(gen, popSize)
 gen.pop3 <- ggplot(gen.pop3, aes(x = gen, y = popSize)) +
          geom_jitter(width = 0.1, height = 0) + labs()
gen.pop3

# population size over generations
```

```{r}
p3dfm2 <- p3mut2 %>%
  select(gen, freqMut)

p3dfm2 <- p3dfm2 %>%   cSplit(p3dfm2, splitCols = "freqMut", sep = '_', direction = "long", fixed = TRUE, drop = FALSE, makeEqual= TRUE, type.convert = TRUE)
# splits freqMut col wherever a '_' occurs, allowing study of frequencies to be made
p3ind <- p3ind %>%   cSplit(p3ind, splitCols = "indvFitness", sep = '_', direction = "long", fixed = TRUE, drop = FALSE, makeEqual= TRUE, type.convert = TRUE)

p3het <- p3het %>%
  cSplit(p3het, splitCols = "Indv_Ho", sep = '_', direction = "long", fixed = TRUE, drop = FALSE, makeEqual= TRUE, type.convert = TRUE) 
p3het<-p3het %>%   mutate(het = (1-p3het$Indv_Ho))
p3het<-p3het %>% select(gen, Indv_Ho, het) %>% 
rename("Hom" = "Indv_Ho")
p3ind <- bind_cols(p3ind, p3het) %>% 
   select(!("gen...11")) %>% 
    rename("gen" = "gen...1")
remove(p3het)

p3dfm2<- p3dfm2 %>% 
  mutate(d.category = case_when(
    freqMut >= 0 & freqMut <= 0.3 ~ "Low",
    freqMut >0.3 & freqMut <= 0.6 ~ "Medium",
    freqMut >0.6 & freqMut < 1.0 ~ "High",
    freqMut == 1.0 ~ "Fixed"
  ))
# categorises m2 frequency into level of effect they may have

p3gen.freqM1 <- ggplot(p3mut1, aes(x = gen, y = Mutcount)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Mutation count (M1)") + theme_classic()

p3gen.freqM2 <- ggplot(p3mut2, aes(x = gen, y = Mutcount)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Deletrious Mutation count (M2)") +theme_classic()

p3gen.freqM3 <- ggplot(p3mut3, aes(x = gen, y = Mutcount)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Beneficial Mutation count (M3)") +theme_classic()

p3gen.freqM1
p3gen.freqM2
p3gen.freqM3
  # plots freq of each mut type. m1 and m3 show common upward trend while m2 is more rollercoastered (any selection would remove del mut whereas m1 and m3 are either neutral or actively beneficial) 
```


```{r}
p3dfm2 <- p3dfm2 %>% 
  mutate(freqMut = as.numeric(freqMut)) %>% 
filter(!is.na(freqMut)) %>% 
       mutate  (indvL = -((-0.01)*(2*0.5* freqMut +(1-(2*0.5)*((freqMut)^2))))
)   
# ^^ works out individual genetic load; 
# below works out additive load vv
# -(sum(del_selCoef*(2*del_domCoef*fMut+(1-(2*del_domCoef)*((fMut)^2)))))
p3dfm2 <- p3dfm2 %>% 
group_by(gen) %>% 
mutate(sum_L_by_gen = (sum(indvL)))

p3m2fixed<- p3dfm2 %>%
  filter(freqMut == 1.0)
#  fixed where freqMut == 1.0,
p3m2seg<- p3dfm2 %>% 
  filter(!(freqMut == 1.0))
#  otherwise segregating.

p3Sum.D.L <- ggplot(p3dfm2, aes(x = gen, y = sum_L_by_gen, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs() +theme_classic()
p3D.L <- ggplot(p3dfm2, aes(x = gen, y = indvL, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs() +theme_classic()
p3D.L
p3Sum.D.L  # if blank then no fixed M2 mutations
# plot of indv genetic load /generation (D.L) and additive load/generation (Sum.D.L)

p3m2Fload.gen <- (ggplot(p3m2fixed, aes(x = gen, y = sum_L_by_gen, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Sum Fixed Load") +theme_classic())

p3m2Sload.gen <- (ggplot(p3m2seg, aes(x = gen, y = sum_L_by_gen, colour = freqMut)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Sum Segregation Load") +theme_classic())

p3m2Fload.gen 
# Drift Load/generation
p3m2Sload.gen
# Segregation Load/generation

```

```{r}
p3tdfm2 <- p3dfm2 %>% 
group_by("d.category") %>% 
count(d.category)
# counts n of each category ("Low", "Medium", and "High")
p3m2.freq <- (ggplot(p3dfm2, aes(x= gen, y= freqMut, colour= d.category)) +
          geom_jitter  (width = 0.1, height = 0) + labs( y = "Mutation frequency") +theme_classic())
p3m2.freq

p3count <- p3dfm2 %>% 
group_by(gen) %>% 
count(d.category)

p3m2.count<- ggplot(p3count, aes(x = gen, y = n, colour = d.category)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Deletrious Mutation count (M2)") +theme_classic()
p3m2.count

# plots the n of categories per gen

```





////////Add sel.coef desc and Combine P1, P2, P3///////////////

```{r}
dfm2 <- dfm2 %>% 
  mutate (sel.coef = 0.014)
p2dfm2 <- p2dfm2 %>% 
  mutate (sel.coef = 0.05)
p3dfm2 <- p3dfm2 %>% 
  mutate (sel.coef = 0.01)

all_dfm2 <- bind_rows(dfm2,  p2dfm2,  p3dfm2)
# put into one table

all_Sum.D.L <- ggplot(all_dfm2, aes(x = gen, y = sum_L_by_gen, colour = sel.coef)) +
          geom_jitter(width = 0.1, height = 0.0) + labs() +theme_classic()
all_Sum.D.L
# ^^ Comparison of Additive Loads

all_m2fixed<- all_dfm2 %>%
  filter(freqMut == 1.0)
#  fixed where freqMut == 1.0,
all_m2seg<- all_dfm2 %>% 
  filter(!(freqMut == 1.0))
#  otherwise segregating.

all_D.L <- ggplot(all_dfm2, aes(x = gen, y = indvL, colour = sel.coef)) +
          geom_jitter(width = 0.1, height = 0) + labs() +theme_classic()
all_D.L
# plot of indv genetic load /generation (D.L)

all_m2Fload.gen <- (ggplot(all_m2fixed, aes(x = gen, y = sum_L_by_gen, colour = sel.coef)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Sum Fixed Load") +theme_classic())

all_m2Sload.gen <- (ggplot(all_m2seg, aes(x = gen, y = sum_L_by_gen, colour = sel.coef)) +
          geom_jitter(width = 0.1, height = 0) + labs(x = "Generation", y = "Sum Segregation Load") +theme_classic())

all_m2Fload.gen   # if blank then no fixed M2 mutations
# Drift Load/generation
all_m2Sload.gen
# Segregation Load/generation
```

```{r}
# all.low.dfm2 <- all_dfm2 %>% 
#   filter(d.category == "Low")


```


```{r}
p1ind <- p1ind %>% 
  mutate (sel.coef = 0.014)
p2ind <- p2ind %>% 
  mutate (sel.coef = 0.05)
p3ind <- p3ind %>% 
  mutate (sel.coef = 0.01)

sel.indv.p1 <- p1ind %>% 
  select(gen, indvFitness, sel.coef)
sel.indv.p2 <- p2ind %>% 
  select(gen, indvFitness, sel.coef)
sel.indv.p3 <- p3ind %>% 
  select(gen, indvFitness, sel.coef)

all_ind <- bind_rows(sel.indv.p1, sel.indv.p2,  sel.indv.p3)
# put into one table


p1.het.sd <- p1ind %>% 
filter(!is.na(het)) %>% 
group_by(gen) %>% 
summarise(
sd_het = sd(het))

p1.het.x <- p1ind %>% 
filter(!is.na(het)) %>% 
group_by(gen) %>% 
summarise(
avg_het = mean(het))

p1.het.stat <- bind_cols(p1.het.x, p1.het.sd) %>% 
select(!("gen...3")) %>% 
  rename("gen" = "gen...1") %>% 
  mutate(sel.coef = 0.014)

p2.het.sd <- p2ind %>% 
filter(!is.na(het)) %>% 
group_by(gen) %>% 
summarise(
sd_het = sd(het))

p2.het.x <- p2ind %>% 
filter(!is.na(het)) %>% 
group_by(gen) %>% 
summarise(
avg_het = mean(het))

p2.het.stat <- bind_cols(p2.het.x, p2.het.sd) %>% 
select(!("gen...3")) %>% 
  rename("gen" = "gen...1") %>% 
  mutate(sel.coef = 0.05)

p3.het.sd <- p3ind %>% 
filter(!is.na(het)) %>% 
group_by(gen) %>% 
summarise(
sd_het = sd(het))

p3.het.x <- p3ind %>% 
filter(!is.na(het)) %>% 
group_by(gen) %>% 
summarise(
avg_het = mean(het))

p3.het.stat <- bind_cols(p3.het.x, p3.het.sd) %>% 
select(!("gen...3")) %>% 
  rename("gen" = "gen...1") %>% 
  mutate(sel.coef = 0.01)



p1.indv.stat <- sel.indv.p1 %>% 
  filter(!is.na(indvFitness)) %>% 
group_by(gen) %>% 
summarise(
sd_fit = sd(indvFitness)) 

p1.indv.stat2 <- sel.indv.p1 %>% 
   filter(!is.na(indvFitness)) %>% 
  group_by(gen) %>% 
  summarise(
    avg_fit = mean(indvFitness)
  )

p1.ind_stat <- bind_cols(p1.indv.stat, p1.indv.stat2) %>% 
select(!("gen...3")) %>% 
  rename("gen" = "gen...1") %>% 
  mutate(sel.coef = 0.014)

p2.indv.stat <- sel.indv.p2 %>% 
  filter(!is.na(indvFitness)) %>% 
group_by(gen) %>% 
summarise(
sd_fit = sd(indvFitness)) 

p2.indv.stat2 <- sel.indv.p2 %>% 
   filter(!is.na(indvFitness)) %>% 
  group_by(gen) %>% 
  summarise(
    avg_fit = mean(indvFitness)
  )

p2.ind_stat <- bind_cols(p2.indv.stat, p2.indv.stat2) %>% 
select(!("gen...3")) %>% 
  rename("gen" = "gen...1") %>% 
  mutate(sel.coef = 0.05)


p3.indv.stat <- sel.indv.p3 %>% 
  filter(!is.na(indvFitness)) %>% 
group_by(gen) %>% 
summarise(
sd_fit = sd(indvFitness)) 

p3.indv.stat2 <- sel.indv.p3 %>% 
   filter(!is.na(indvFitness)) %>% 
  group_by(gen) %>% 
  summarise(
    avg_fit = mean(indvFitness)
  )

p3.ind_stat <- bind_cols(p3.indv.stat, p3.indv.stat2) %>% 
select(!("gen...3")) %>% 
  rename("gen" = "gen...1") %>% 
  mutate(sel.coef = 0.01)
```


```{r}
all.ind.stat <- bind_rows(p1.ind_stat, p2.ind_stat, p3.ind_stat)

p1.w0 <- all.ind.stat %>% 
  filter(gen == 2) %>% 
  filter(sel.coef == 0.014)
p2.w0 <- all.ind.stat %>% 
  filter(gen == 2) %>% 
  filter(sel.coef == 0.05)
p3.w0 <- all.ind.stat %>% 
  filter(gen == 2) %>% 
  filter (sel.coef == 0.01)

p1.ind_stat <- p1.ind_stat %>% 
 group_by(sel.coef) %>% 
 mutate (wt.w0 = (p1.ind_stat$avg_fit/(p1.w0$avg_fit)))
p2.ind_stat <- p2.ind_stat %>% 
 group_by(sel.coef) %>% 
 mutate (wt.w0 = (p2.ind_stat$avg_fit/(p2.w0$avg_fit)))
p3.ind_stat <- p3.ind_stat %>% 
 group_by(sel.coef) %>% 
 mutate (wt.w0 = (p3.ind_stat$avg_fit/(p3.w0$avg_fit)))
# mutate (wt.w0 = (avg_fit/(avg.fit (gen row1))))

all.ind.stat <- bind_rows(p1.ind_stat, p2.ind_stat, p3.ind_stat)

# d.W = Wt/W0 = change in fitness (y)// generation (x)
d.W_gen <- (ggplot(all.ind.stat, aes(x = gen, y = wt.w0, colour = sel.coef)) + geom_jitter  (width = 0.1, height = 0) + labs( x = "Generation", y = "Wt/W0") +theme_classic())
d.W_gen
```

```{r}
all.het.stat <- bind_rows(p1.het.stat, p2.het.stat, p3.het.stat)

p1.h0 <- all.het.stat %>% 
  filter(gen == 2) %>% 
  filter(sel.coef == 0.014)
p2.h0 <- all.het.stat %>% 
  filter(gen == 2) %>% 
  filter(sel.coef == 0.05)
p3.h0 <- all.het.stat %>% 
  filter(gen == 2) %>% 
  filter (sel.coef == 0.01)

p1.het.stat <- p1.het.stat %>% 
 mutate (ht.h0 = (p1.het.stat$avg_het/(p1.h0$avg_het)))
p2.het.stat <- p2.het.stat %>% 
 mutate (ht.h0 = (p2.het.stat$avg_het/(p2.h0$avg_het)))
p3.het.stat <- p3.het.stat %>% 
 mutate (ht.h0 = (p3.het.stat$avg_het/(p3.h0$avg_het)))
# mutate (wt.w0 = (avg_fit/(avg.fit (gen row1))))

all.het.stat <- bind_rows(p1.het.stat, p2.het.stat, p3.het.stat)

# d.W = Wt/W0 = change in fitness (y)// generation (x)
d.H_gen <- (ggplot(all.het.stat, aes(x = gen, y = ht.h0, colour = sel.coef)) +  geom_jitter  (width = 0.1, height = 0) + labs( x = "Generation", y = "Ht/H0") +theme_classic())
d.H_gen

```
```{r}

all.ind.stat <- bind_cols(all.ind.stat, all.het.stat$ht.h0) %>% 
  rename(ht.h0 = ...6)
d.W_d.H <- (ggplot(all.ind.stat, aes(x = ht.h0, y = wt.w0, colour = sel.coef)) +
         geom_jitter  (width = 0.1, height = 0) + labs( x = "Ht/H0", y = "Wt/W0") +theme_classic())
d.W_d.H
```

  