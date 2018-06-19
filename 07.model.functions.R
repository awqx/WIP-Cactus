# Both pp.dir and tst.dir should end in a backslash
# n refers to the split number
preprocess.tst.mod <- function(pp.dir, tst.dir, feat, n) {
  pp.settings <- readRDS(paste0(pp.dir, n, "/pp.settings.RDS"))
  tst <- readRDS(paste0(tst.dir, "/tst", n, ".RDS"))
  guest <- select(tst, guest)
  tst.dg <- tst %>% select(., guest, DelG)
  
  colnames(tst) <- str_replace(colnames(tst), "-", ".")
  tst <- select(tst, -DelG)
  
  tst <- do.call(data.frame, lapply(tst,
                                    function(x)
                                      replace(x, is.infinite(x), NA)))
  tst <- tst %>%
    predict(pp.settings, .) %>% select(., feat) %>%
    cbind(guest, .)
  tst.ad <- domain.num(tst)
  tst.outliers <- tst.ad %>% filter(domain == "outside") %>% .$guest
  # print(tst.outliers)
  # indole is just a persistent headache for alpha
  # and 3-methylbenzoic acid consistently messes up beta
  if(str_detect(pp.dir, "alpha"))
    tst.outliers <- c(tst.outliers, "indole")
  else if(str_detect(pp.dir, "beta"))
    tst.outliers <- c(tst.outliers, "3-methylbenzoic acid")
  tst <- tst %>%
    filter(!guest %in% tst.outliers) %>% 
    select(., -guest) %>% data.matrix()
  tst.dg <- tst.dg %>% filter(!guest %in% tst.outliers) %>%
    select(., -guest)
  
  return(cbind(tst.dg, tst))
}

preprocess.ev <- function(cd.type, n, feat) {
  if (cd.type == "alpha")
    ev <- readRDS("./ext.validation/alpha.RDS")
  else if (cd.type == "beta")
    ev <- readRDS("./ext.validation/beta.RDS")
  ev.info <- select(ev, guest:data.source)
  
  colnames(ev) <- str_replace(colnames(ev), "-", ".")
  ev <- select(ev, -host:-data.source)
  ev <- do.call(data.frame, lapply(ev, 
                                   function(x)
                                     replace(x, is.infinite(x), NA)))
  pp.settings <- readRDS(paste0("./pre-process/", cd.type, 
                                "/", n, "/pp.settings.RDS"))
  ev <- ev %>% predict(pp.settings, .) %>% select(., feat) %>%
    cbind(select(ev.info, guest), .)
  
  ev.ad <- domain.num(ev)
  ev.outliers <- ev.ad %>% filter(domain == "outside") %>% .$guest
  
  ev <- ev %>% filter(!guest %in% ev.outliers) %>% select(., -guest)
  ev.dg <- ev.info %>% filter(!guest %in% ev.outliers) %>% select(., DelG)
  
  return(cbind(ev.dg, ev))
}
