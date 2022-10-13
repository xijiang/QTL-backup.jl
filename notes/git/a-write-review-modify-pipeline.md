---
title: A blog-review-publish pipeline
author: Xijiang Yu
institute: Norwegian University of Life Sciences
date: Oct. 5, 2022
theme: Warsaw
---

## Objectives
- To demonstrate MarkDown writing of a blog
  - can also be a paper
- To demonstrate to keep a history of this writing
- To demonstrate cooperation (review) of the same repo

## The beginning
- Sam want to write a quality blog
- To make it of better quality
- He is asking Tom to review
- Sam and Tom share a same server `bridge` with some arbitrary IP
- They decided to use git and markdown tools for this blog

## Preparation
```bash
# On server bridge
cd /mnt/share/md  # a shared directory that both Sam and Tom has r/w rights
mkdir blog
cd blog
git init --bare   # this directory will server solely for a repo, no working space
```

## Sam initiates the blog
- Clone repo from the `bridge`
```bash
# in WSL of Sam's local laptop
mkdir work
cd work
git clone sam@bridge:/mnt/share/md  # a md directory appeared now in dir work
cd md
touch myblog.md
git add myblog.md
```

## Sam pop in more contents
- Outlining of the blog
```markdown
---
title: A blog
author: Sam
date: \today
---
## Abstract
## Introduction
## Materials and methods
## Results
## Referencess
```

## Sam is satisfied for the structure
```bash
git add myblog.md  # again, to stage it and ready for commit
git commit
```

- An editor is activated asking Sam to write message for this commit
- Sam wrote a title line, an empty line, and a paragraph, which is a quite good habbit

| |
|--|
|Structure of the blog|
||
|The blog is planning big|

- Sam then
```bash
# push the local modifications to bridge for Tom to review
git push origin main
```

## Tom starts to review
```bash
# In WSL of Tom's laptop
mkdir work
cd work
git clone tom@bridge:/mnt/share/md
```

- Tom then edit the `myblog.md` in `md` directory
- He found there should be an *conclusion* section
- There is also a typo, an extra *s* in the last line

## Tom made the correction, commit and push
```bash
# after modification
git add myblog.md # ready to be committed again
git commit        # message is also a title line and some explanation
git push origin main
```

- Then Tom notified Sam to check

## Sam works on the modified repo
```bash
# in his `~/work/md/`
git pull  # all the modifications by Tom is now local to Sam
git log   # now Sam can see to snapshots (commits) in history
git diff HEAD HEAD~1 # to see what was modified.
```

- This pipeline can be repeated between Sam and Tom until they a satisfied with the blog.
- Finally, Sam and Tom realize this blog can be published in some journal

```bash
pandoc myblog.md -o thesis.docx  # or
pandoc myblog.md -o thesis.pdf
```

## Other issues
- It is better to write the blog one sentenc or shorter a line.
- When Sam and Tom edited the same file and both want to push
- The later push will be rejected, because of confilictions
  - a 3-way merge is required, which is OK
  - different people work on different sections would be better, and then

```bash
pandoc -o myblog.pdf \
	intro.md \
	method.md
```

- The remote server can be https://github.com
  - a fork $\to$ pull-request $\to$ code-reviews $\to$ merge pipe line is similar
