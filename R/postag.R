#' do part-of-speech tagging
#' 
#' Tag text(s) with part-of-speech tagger
#' @export
postag <- function(x, ...) {
  UseMethod("postag")
}

#' @rdname postag
#' @export
postag.character <- function(x, model = NULL, method = c("perceptron", "hmm")) {
  if (method == "perceptron") {
    stop("Perceptron POS tagger is not implemented yet")
  }
  else if (method == "hmm") {
    if (is.null(model)) {
      # require(quantedaData)
      # model <- data(hmmModel, package = "quantedaData")
    }
    viterbi(x, transition = model[[2]], observation = model[[1]])
  }
  else {
    stop(method, " not implemented in postag().")
  }
}

#' @rdname postag
#' @export
postag.corpus <- function(x, method = c("perceptron", "hmm")) {
  postag(texts(x))
}

#' train hidden Markov model(HMM)
#' 
#' Train hidden Markov model(HMM) which maximises the joint probability of
#' token and tag sequence. This is done by aggregating the frequencies of tags
#' being assigned to a specific token and the frequencies of transitions between
#' tags. The function takes labelled list of sentences as an input, where each
#' element of a sentence is a list consisting of token and tag.
#' @param labeled a list of sentences where each element of a sentence is a list
#' consisting of token and tag
#' @param alpha adjustable parameter used for smoothing (alpha must be between 0 
#' and 1)
#' @return A list of two data.tables with token observation and tag transition 
#' probabilities
#' @examples 
#' @export
trainHMM <- function(labeled, alpha = 0.1) {
  if (is.list(labeled) && is.null(dim(labeled))) {
    tokens <- unlist(lapply(labeled, "[[", 1))
    tags <- unlist(lapply(labeled, "[[", 2))
  }
  else {
    stop("Labeled sequence must be a list.")
  }
  if (alpha < 0 | alpha > 1) {
    stop("Alpha must be between 0 and 1.")
  }
  tokentag <- mleTokenTag(tokens, tags, alpha)
  dt <- data.table()
  for (i in 1:length(labeled)) {
    tmp <- data.table("tag" = unlist(lapply(labeled[i], "[[", 2), use.names = FALSE))
    tmp <- tmp[, previous := shift(tag)]
    dt <- rbind(dt, tmp)
  }
  tagtrans <- mleTagTransition(dt)
  model  <- list(tokentag, tagtrans)
  return(model)
}

#' tag text using Viterbi algorithm
#' 
#' Tag text using hidden Markov model(HMM) with Viterbi algorithm for decoding. 
#' Use either pre-trained models from quantedaData package or train your own model
#' with \code{train()} function.
#' @param txt character string
#' @param transition a table of tag transition probabilities created with 
#' \code{train()} function.
#' @param  observation a table of state observation likelihoods created with
#' \code{train()} funcion
#' @return A vector of part-of-speech tags
#' @references Jurafsky, D. & Martin J. (2009). Speech and Language Processing.
#' @examples 
#' @export
viterbi <- function(txt, transition = NULL, observation = NULL) {
  if (is.null(transition)) {
    stop("Transition probability matrix must be specified.")
  }
  if (is.null(observation)) {
    stop("State observation likelihood matrix must be specified.")
  }
  
  tokens <- unlist(tokenize(txt))
  
  states <- as.character(unique(transition[,tag]))
  
  # Viterbi path probability trellis
  pathmat <- data.table(matrix(0, nrow = length(states), ncol = length(tokens)))
  
  # Initialization
  for (s in 1:length(states)) {
    transprob <- transition[tag == states[s] & previous == "NA", prob]
    statobsprob <- getStateObsProb(tokens[1], states[s], observation)
    pathmat[s, V1 := transprob * statobsprob]
  }
  # Recursion
  for (t in 2:length(tokens)) {
    for (s in 1:length(states)) {
      viterbiprob <- max(pathmat[,t-1] * transition[tag == states[s], prob])
      statobsprob <- getStateObsProb(tokens[t], states[s], observation)
      pathmat[s, paste0("V", as.character(t)) := viterbiprob * statobsprob]
    }
  }
  # Termination
  pos <- states[apply(pathmat, 2, which.max)]
  pos
}

# Maximum Likelihood estimator for calculating the probability  of observing 
# token given tag (emission probability). Adjustable parameter alpha is used
# to implement simple Laplace(Lidstone) smoothing for out-of-vocabulary(OOV) tokens. 
mleTokenTag <- function(tokens, tags, alpha = 0.1) {
  if (typeof(tokens) != "character") {
    stop("Vector of tokens must be a character vector.")
  }
  if (typeof(tags) != "character") {
    stop("Vector of tags must be a character vector.")
  }
  if (length(tokens) != length(tags)) {
    stop("Length of tokens vector isn't equal to the length of tags vector.")
  }
  uniquetags <- unique(tags)
  dt <- data.table("token" = tokens, "tag" = tags)
  dt <- dt[, 
           .(countTokenTag = .N), 
           by = list(token, tag)
           ]
  temp <- data.table("token" = rep("OOV", length(uniquetags)),
                     "tag" = uniquetags,
                     "countTokenTag" = 0)
  dt <- rbindlist(list(dt, temp))
  dt <- dt[,
           c("totalTag", "uniqueWords") := list(sum(countTokenTag), uniqueN(.SD)),
           by = tag
           ][,
              prob := (countTokenTag + alpha)/(totalTag + alpha * uniqueWords)
            ]
  return(dt[, .(token, tag, prob)])
}

# Retrieve probability of observing token given tag (state observation probability).
# Currently implemented smoothing method (Laplace) doesn't take into account the
# morphology of unknown words or differing prior probabilities for different 
# parts of speech.
getStateObsProb <- function(state, obs, tab) {
  prob <- tab[token == state & tag == obs, prob]
  if (length(prob) == 0) {
    prob <- 0
    prob <- tab[token == "OOV" & tag == obs, prob]
  }
  return(prob)
}

# Maximum Likelihood estimator for calculating the between-tag transition.
# The starting point is coded as NA, which gives the probability of beginning
# a sentence with a given tag.
mleTagTransition <- function(dttags) {
  dt <- dttags[,
               totalPreviousTag := .N,
               by = previous
               ][,
                 .(countTransition = .N),
                 by = list(tag, previous, totalPreviousTag)
                 ][,
                   prob := countTransition/totalPreviousTag
                   ][is.na(previous),
                     previous := "NA"
                     ]
  dt <- buildTransitionProbabilityMatrix(dt)
  return(dt[,.(tag, previous, prob)])
}

# Populate transition probability matrix with unseen tag transitions.
buildTransitionProbabilityMatrix <- function(dt) {
  # Generate all permutations of tags
  uniquetags <- unique(dt[, previous])
  transitions <- data.table(expand.grid(tag = uniquetags, previous = uniquetags))
  setkey(transitions, tag, previous)
  setkey(dt, tag, previous)
  
  seen <- unique(transitions[dt[,.(tag,previous)], which = TRUE])
  transitions[seen, prob := dt[,prob]]
  transitions[!seen, prob := 0]
  return(transitions)
}
