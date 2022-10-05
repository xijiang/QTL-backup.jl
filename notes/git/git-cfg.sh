#!/usr/bin/env bash
git config --global user.name "Xijiang Yu"
git config --global user.email xijiang@users.noreply.github.com
git config --global core.editor emacs # if you like Emacs more, default vi
git config --global merge.tool meld
git config --global alias.ss 'status -s'
git config --global alias.lg "log --color --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)<%an>%Creset' --abbrev-commit --branches"
git config --global alias.tg "log --date-order --graph --tags --simplify-by-decoration --pretty=format:'%ai %h %d'"
git config --global init.defaultBranch main
