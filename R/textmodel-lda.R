#' Convert a \link{dfm} into the format needed by \pkg{lda}
#' 
#' Convert a quanteda \link{dfm} object into the indexed format required by the
#' topic modelling package \pkg{lda}.
#' @param d A \code{\link{dfm}} object
#' @return A list with components "documents" and "vocab" as needed by 
#'   \link[lda]{lda.collapsed.gibbs.sampler}
#' @import topicmodels
#' @export
#' @examples
#' mycorpus <- subset(inaugCorpus, Year>1970)
#' d <- dfm(mycorpus, stopwords=TRUE)
#' d <- trim(d, minCount=5, minDoc=3)
#' td <- dfm2ldaformat(d)
#' if (require(lda)) {
#'     tmodel.lda <- lda.collapsed.gibbs.sampler(documents=td$documents, 
#'                                               K=10,  
#'                                               vocab=td$vocab,
#'                                               num.iterations=50, alpha=0.1, eta=0.1) 
#'     top.topic.words(tmodel.lda$topics, 10, by.score=TRUE) # top five words in each topic
#' }
dfm2ldaformat <- function(d) {
    tmDTM <- dfm2tmformat(d)
    return(topicmodels::dtm2ldaformat(tmDTM))
}

#' latent Dirichlet allocation text model
#' 
#' \code{textmodel_lda} estimates the parameters of Blei et. al. (2003)
#' This function is a wrapper to \link[topicmodels]{LDA} and
#' \link[topicmodels]{CTM} in \pkg{topicmodels}.
#' @param x the dfm on which the model will be fit
#' @param model the class of lda-type model to fit.  \code{lda} for classic 
#'   latent Dirichlet allocation; \code{ctm} for Blei and Rafferty's (2007) 
#'   correlated topic model.
#' @param k number of topics  
#' @param smooth a smoothing parameter for word counts, defaults to 0
#' @param meta metadata passed to \link[stm]{stm} if fitting a structural topic
#'   model
#' @param ... additional arguments passed to \link[topicmodels]{LDA}
#' @references Blei, David M., Andrew Y. Ng, and Michael I. Jordan. 2003. "Latent 
#'   Dirichlet Allocation." \emph{The Journal of Machine Learning Research} 3: 
#'   993-1022.
#'   
#'   Blei, David M, and John D. Lafferty. 2007. "A Correlated Topic Model of 
#'   Science." \emph{The Annals of Applied Statistics} 1(1): 17-35.
#'   
#'   Roberts, M., Stewart, B., Tingley, D., and Airoldi, E. (2013) "The
#'   structural topic model and applied social science." In \emph{Advances in Neural
#'   Information Processing Systems Workshop on Topic Models: Computation,
#'   Application, and Evaluation}. \url{http://goo.gl/uHkXAQ}
#' @import topicmodels stm
#' @author Kenneth Benoit
#' @examples 
#' \dontrun{
#' data(sotuCorp, package="quantedaData")
#' SOTUCorpus <- sotuCorp
#' presDfm <- dfm(subset(SOTUCorpus, year>1960), stopwords=TRUE, stem=TRUE)
#' presDfm <- trim(presDfm, minCount=5, minDoc=3)
#' presLDA <- textmodel_lda(presDfm, k=10)
#' #require(topicmodels)  # need this to access methods below
#' terms(presLDA, k=10)  # top 10 terms in each topic
#' topics(presLDA)       # dominant topics for each document
#' presCTM <- textmodel_lda(presDfm, model="ctm", k=10)
#' terms(presCTM, k=10)  # top 10 terms in each topic
#' topics(presCTM)       # dominant topics for each document
#' 
#' # fit a structural topic model
#' if (require(stm)) {
#'   gadarianCorpus <- corpus(gadarian$open.ended.response, docvars=gadarian[, 1:3],
#'                            source="From stm package, from Gadarian and Albertson (forthcoming)")
#'     gadarianDfm <- dfm(gadarianCorpus, stopwords=TRUE, stem=TRUE)
#'     gadarianSTM <- textmodel_lda(gadarianDfm, "stm", k=3, 
#'                                  prevalence = ~treatment + s(pid_rep), 
#'                                  data = docvars(gadarianCorpus))
#'     summary(gadarianSTM)
#' }
#' }
#' @import topicmodels stm
#' @export
textmodel_lda <- function(x, model=c("lda", "ctm", "stm"), k, smooth=0, meta=NULL, ...) {
    model <- match.arg(model)
    if (!is.dfm(x))
        stop("supplied data must be a dfm object.")
    if (smooth) 
        x <- x + smooth  # smooth by the specified amount if > 0
    if (model=="lda") {
        dataSTM <- dfm2tmformat(x)
        #require(topicmodels)
        fittedlda <- topicmodels::LDA(dataSTM, k=k, ...) 
    } else if (model=="ctm") {
        dataSTM <- dfm2tmformat(x)
        fittedlda <- topicmodels::CTM(dataSTM, k=k, ...)   
    } else if (model=="stm") {
        stmdata <- dfm2stmformat(x)
        fittedlda <- stm::stm(stmdata$documents, stmdata$vocab, K=k, verbose=FALSE, ...)
    }
    return(fittedlda)
}


# rdname predict.textmodel
# param ... additional arguments passed to \link[austin]{predict.wordfish}
# references LBG (2003); Martin and Vanberg (2007)
# return \code{predict.wordfish} is a wrapper to \link[austin]{predict.wordfish}, and has the same returns, 
# # except that the object is also classes as "wordfish".
# author Ken Benoit, wrapper around Will Lowe's \pkg{austin} code.
# export
# predict.wordfish <- function(object, newdata=NULL, ...) {
#     
#     result <- austin::predict.wordfish(object, newdata, ...)
#     class(result) <- c("wordfish", "data.frame")
#     return(result)
# }



