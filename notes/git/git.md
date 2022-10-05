---
title: Best practices of GIT
author: Xijiang Yu
institute: Nowegien University of Life Sciences
date: Oct. 5, 2022
theme: Warsaw
---

## What is not covered
- How to setup an account on [GitHub](https://git.com).
  - You should
  - Better with private-public key pairs
  - Alternatives: [BitBucket](https://bitbucket.org/dashboard/overview), and [GitLab](https://about.gitlab.com/).
- Why git is the best version control system (VCS)
- The mechenisms behind `git`

## Bonus if we have time
- Create ssh key pairs
- Using `tmux`

## Installation
- Windows
  - [64 bit version](https://github.com/git-for-windows/git/releases/download/v2.38.0.windows.1/Git-2.38.0-64-bit.exe)
  - WSL (recommended)
```sh
# in windows powershell
wsl --install
wsl --list --online
wsl --install -d Ubuntu # Install Ubuntu 20.04 LTS
# maybe after reboot and then in WSL -> Ubuntu
sudo apt update
sudo apt install git meld # emacs etc...
```
- Linux
  - You will have an account on an intranet server
  - Everything is ready there.
  
## Configuration
```bash
# run `git-cfg.sh` together with these slides.
git config --global user.name "Xijiang Yu"
git config --global user.email xijiang@users.noreply.github.com
git config --global core.editor emacs # if you like Emacs more, default vi
git config --global merge.tool meld
git config --global alias.ss 'status -s'
git config --global alias.lg "log --color --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)<%an>%Creset' --abbrev-commit --branches"
git config --global alias.tg "log --date-order --graph --tags --simplify-by-decoration --pretty=format:'%ai %h %d'"
git config --global init.defaultBranch main
```

## The structure of a `git` directory
```bash
# In one of your working directories
mkdir init
cd init
git init
tree .git

## Alternatively
mkdir bare
cd bare
git init --bare
tree .
# `bare` will only serve as a repo
```

## The many states of a file in a git system.
```bash
mkdir -p ~/test/git
cd ~/test/git
git init     # starting a repo with working space
touch me.md  # created an empty untracted file
git add me.md # this file is now tracked and staged
# Add some contents into the file with any text editor, or
echo I am John. I am learning git >> me.md # now the file is modified
git status
git add me.md # staged again
# then more modification
git commit -am 'my first version controlled file' # commited
git log
```

## Branching and collaboration
- Create a repo that have introducitons of every party in the club
- Also practice some MarkDown

```bash
# Preparation on the intranet server
cd /home/share
mkdir git-test
cd git-test
git init  # and then some files by me

# In your working directory
git clone /home/share/git-test .
git checkout -b your-name
touch your-name.md
# make modifications and finally commit
git push # send modifications (as if) to the server side

# See the merge precess
```

## What's next
- Create a GitHub account
- Create a (private) repo to save all you useful scripts and notes
- git clone a local copy of this repo
- Add contents, modifications and snapshots
- Push to github
- Stick to the repo to see if you can still use it after 10 years.
- More will be talked about in the Julia section.

## Further reading
- Free courses
- https://git-scm.com/book/en/v2
- And many more.
