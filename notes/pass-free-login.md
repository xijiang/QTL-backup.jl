---
title: Instructions for password free login to remote servers
author: Xijiang Yu
institute: Norwegian University of Life Sciences
date: Oct. 5, 2022
---
## The problem
- Password login to a remote server is tedious and not safe
  - [A 8 mixed case + numbers password cost only $155 to decode in 2020](https://turingpoint.de/en/blog/how-much-does-it-cost-to-crack-your-password/)
- Key pairs are safer and hassle free

## Generate a key pair
- In Widnows 10/11 WSL
- In a Linux/MacOS terminal, run the following commands
```bash
sudo apt install openssh # for WSL if necessary
ssh-keygen
# the key pairs will be saved in ~/.ssh by default
# ~/.ssh must of property 700
chmod 700 ~/.ssh # if necessary
# then input a passphrase, which is a phrase, 
# it can be much more complicated than a (pass)word
ssh-copy-id user@remote # to copy the public key to the remote server
ssh user@remote # now you're asked for the phrase, instead of word to login
```

## Private key management
- You will have to enter the more complicated phrase to login to remote now
- On a Linux system, when you logged in locally via GUI, e.g., gnome
  - system will ask the passphrase for the private key the first time you use it
  - it also ask if you want to store the phrase
  - you will then never be asked for the phrase and password for remotes again if using GUI local
- for WSL, and remotes terminal
```bash
sudo apt install keychain
# insert below line into ~/.bashrc
eval $(keychain --eval ~/.ssh/id_rsa)
```
- Then you will be asked for phrase once for the session, until it ends
- If `keychain` is not available, for a pass-free session, run
```bash
eval $(ssh-agent)
ssh-add # to let the agent manage the phrase for this session
```

## Use the key pair for github
- Register a GitHub account
- Navigate to https://github.com/settings/ssh/new
  - copy `~/.ssh/id_rsa.pub` contents and paste them into the `key` input box
  - give it a proper title also
- You can then use command lines for your codes backup.

