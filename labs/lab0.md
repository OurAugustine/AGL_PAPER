# 0. Setting the Stage
It is recommended to install Python in your local machine. You can download the latest version of Python from [Anaconda](https://www.anaconda.com/). Please read the installation documents on the Anaconda website for your own Operation System (i.e. Windows OS, Mac OS, or linux).

## Discord
We will use [Discord](https://discord.com/) as our main communication channel. Please join the [Discord server](https://discord.gg/VaaJHznX) for this course. We will use the `#labs` channel for all labs.
Please use this adventage to get to know your classmates and ask questions about the course.

**NOTE**: The Discord server is not set up yet. I will update the link here once it is set up.

## ChatGPT and Bard
ChatGPT and Bard are two chatbots that can help you with your work.
Please use the advantage of these chatbots when you have questions at the first place.
They are the most accessible tools than any other human beings. 

### Prompt engineering
Here are strategies for getting better results from chatGPT, which is gathered from the [OpenAI docs](**Reference**: https://platform.openai.com/docs/guides/prompt-engineering/six-strategies-for-getting-better-results). Other strategies require GPT4 models, which is under a subscription.

#### 1. Write clear instructurions
As the model cannot read your mind. You need to write clear instructions for the model to understand what you want.
If outputs are too long, ask for brief replies. If outputs are too simple, ask for more details and examples.
If you don't like the format of the output, provide a sample output you'd like to see.
The less the model has to guess, the better the results will be.

**Try**: Write a clear question for ChatGPT to help you install Python according to your OS.

#### 2. Split complex tasks into simpler subtasks
Just as it is good practice in software engineering to decompose a complex system into a set of modular components,
the same is true of tasks submitted to a language model. Complex tasks tend to have higher error rates than simpler tasks.
Further more, complex tasks can often be re-defined as a workflow of simpler tasks in which the outputs of earlier tasks are used to construct the inputs to later tasks.

#### 3. Give the model time to "think"
Models make more reasoning errors when trying to answer right away, rather than taking time to work out an answer. 
Asking for a "chain of thought" before an answer can help the model reason its way toward correct answers more reliably.

#### 4. Test changes systematically
Improving performance is easier if you can measure it. In some cases a modification to a prompt will achieve better performance on a few isolated examples but lead to worse overall performance on a more representative set of examples. Therefore to be sure that a change is net positive to performance it may be necessary to define a comprehensive test suite (also known an as an "eval").

## Git and Gitlab
We will use [Git](https://git-scm.com/) and [Gitlab](https://git.uark.edu/) for submitting your assignments.
The submission format of your assignments will be text and python files. 
[Github](https://github.com) is another popular version control repository. However, we will use Gitlab for this course.

**Note:** You need to install git in your local machine. You can download the latest version of git from [here](https://git-scm.com/downloads).

### GitLab SSH Key
You need to generate a ssh key for your Gitlab account. You can follow the instructions [here](https://docs.gitlab.com/ee/user/ssh.html) to generate a ssh key for your Gitlab account. You can use the following command to generate a ssh key.

#### Generating an SSH Key

1. **On Windows:**
   - Open Git Bash (you can download it from [Git for Windows](https://gitforwindows.org/) if you don't have it).
   - Use the command:
     ```bash
     ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
     ```
   - When prompted to "Enter a file in which to save the key," press Enter to accept the default file location.
   - Enter a secure passphrase when prompted.

2. **On MacOS:**
   - Open the Terminal.
   - Run the same command:
     ```bash
     ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
     ```
   - Follow the prompts as above to choose a save location and passphrase.

#### Adding the SSH Key to GitLab

1. **Add the Key to GitLab:**
   - Go to GitLab and sign in to your account.
   - Click on your profile picture in the top right corner, then click on "Settings".
   - In the User Settings sidebar, click on "SSH Keys".
   - Paste your key into the "Key" field.
   - Add a title in the "Title" field to help you identify the key (e.g., "MacBook Pro" or "Windows Desktop").
   - Click "Add key".

Once you've added your SSH key to GitLab, you should be able to perform Git operations over SSH, which is more secure and often more convenient than using HTTPS, especially for frequent Git users.

Remember, the email you use in the `ssh-keygen` command should be the same one associated with your GitLab account. This helps in easily identifying the key's purpose and also in troubleshooting, should any issues arise with your GitLab connection.

### **Task 1: Create a Gitlab account and a repository for this course.**
- **Step 1**: Use your UARK email and password to login to [Gitlab](https://git.uark.edu/).
- **Step 2**: Go to the website [MATH-499V-599V-MFML](https://git.uark.edu/jiahuic/math-499v-599v-mfml) and click the button `Fork` on the top right corner.
- **Step 3**: Name your forked repository as `MFML-{First Name}-{Last Name}`. For example, `MFML-John-Doe`.
- **Step 4**: Clone Your Forked Repository. You can use the following command to clone your forked repository to your local machine. You can try the following code in your Powershell (WinOS) or terminal (MacOS) or VSCode (both).
```bash
git clone https://git.uark.edu/{your-username}/MFML-{First Name}-{Last Name}.git
```
- **Step 5**: Add me as a **Developer** to your forked repository. You can do this on your repository page. Click `Manage` on the left side bar and choose `Members`. Click `Invite members` on the rightup corner of the webpage. In the search panel, check `jiahuic` or `jiahuic@uark.edu` and select the role as `Developer` and clike `invite`.
<!-- - **Step 6**: -->
<!-- - **Step 7**: -->

### **Task 2: Practice Git commands. `git add`, `git commit`, `git push`, `git pull`**
The second task 2 is to practice the basic git commands. You can use the following commands to practice git commands. At the beginning of our course, it is enough to use the following Git commands. As we go deeper of this course, we will revisit more git skill.

- **Step 1**: Create a new folder at `./labs/` and name it as `lab0`. You can use the following command to create a new folder.
```bash
mkdir ./labs/lab0
```
- **Step 2**: Create a markdown file at `./labs/lab0/` and name it as `README.md`. You can use the editor of your choice to create this file. Or, you can try the following commend to create this file.
```bash
touch ./labs/lab0/README.md
```
Or use your editor to create this file.
- **Step 3**: Editor the `README.md` file. On the first line, type `# Lab 0`. On the second line, type `This is a lab for practicing Git commands.`. On the third line, type `Please type your name below.`. On the fourth line, type your name. For example, `John Doe`.
- **Step 4**: Add the `README.md` file to the staging area. You can use the following command to add the `README.md` file to the staging area.
```bash
git add ./labs/lab0/README.md
```
- **Step 5**: Commit the changes. You can use the following command to commit the changes.
```bash
git commit -m "Add README.md file"
```
- **Step 6**: Push the changes to your forked repository. You can use the following command to push the changes to your forked repository.
```bash
git push origin main
```
- **Step 7**: Check the changes on your forked repository at the GitLab website. You can see the changes on your forked repository at the GitLab website.

Overall, the steps follow that you make changes to your local repo, add the changes to the staging area, commit the changes, and push the changes to your forked repository. You can check the changes on your forked repository at the GitLab website.

NOTE: You can have the full credits after you finish the first two tasks.

### **Task 3: Practice Fork repository and Pull request**
The third task is to practice the fork repository and pull request. You can use the following commands to practice fork repository and pull request as I will update the code in the center repository.
This case happens when your forked repository is both behind and ahead of the center repository, it means you have made changes that aren't in the upstream repository (ahead) and there are changes in the upstream repository that you don't have (behind). To synchronize your repository with the upstream while keeping your own commits, you can use a combination of `git fetch`, `git merge`, or `git rebase`.
To check if the central (upstream) repository has been updated since you forked or last updated your repository, you can follow these steps:

1. **Configure a Remote for the Upstream Repository:**
   If you haven't already configured a remote for the upstream repository (the original one that you forked), you'll need to do that first. You can do this by running:
   ```bash
   git remote add upstream URL_OF_CENTRAL_REPO
   ```
   Replace `URL_OF_CENTRAL_REPO` with the actual URL of the central repository. `URL_OF_CENTRAL_REPO = https://git.uark.edu/jiahuic/math-499v-599v-mfml.git` for this course.

2. **Fetch the Latest Changes from Upstream:**
   Now, fetch the changes from the upstream repository. This command won't merge any changes into your local repository; it will just fetch the updates for review. Pull updates from the upstream repository.
   ```bash
   git fetch upstream
   git merge upstream/main
   ```
3. **Check the Status Against Upstream:**
   To see if your local branch is behind the upstream repository's branch, you can use:
   ```bash
   git status
   ```
   Or, to see more detailed information, compare your branch with the upstream branch:
   ```bash
   git diff --stat HEAD upstream/main
   ```
   Replace `main` with whatever branch you are tracking (e.g., `master` if the main branch is named `master`). We use `main` for this course.
   - Using `git status` after merging helps to ensure that your local branch is up-to-date with the upstream branch.
   - `git diff --stat HEAD upstream/main` is useful for seeing a summary of what changes have been incorporated from the upstream.

4. **Viewing the Commit Log:**
   If you want to see the actual commits that are in the upstream repository but not in your repository, you can use:
   ```bash
   git log HEAD..upstream/main --oneline
   ```
   This command will list the commits present in the upstream `main` branch but not in your current branch. Again, replace `main` with the appropriate branch name if different.

5. **Merge the Upstream Changes into Your Repository**
   Now, you have a couple of choices on how to integrate these changes: [merging](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging) or [rebasing](https://git-scm.com/docs/git-rebase).

   - Option 1: Merging
       Merging is a non-destructive operation and will keep the history of your changes and the upstream changes. 
       
       First, switch to your main branch (if you're not already on it):
       ```bash
       git checkout main
       ```
       
       Then, merge the changes from the upstream repository:
       ```bash
       git merge upstream/main
       ```
       
       This command will create a new commit in your history that merges the two histories. If there are any merge conflicts, you'll need to resolve them. Use `git diff` in step 3 to check the difference.

   - Option 2: Rebasing
       Rebasing will place your changes on top of the upstream changes, resulting in a cleaner, linear history.
       
       First, switch to your main branch:
       ```bash
       git checkout main
       ```
       
       Then, rebase your changes onto the upstream branch:
       ```bash
       git rebase upstream/main
       ```
       
       During a rebase, you might encounter conflicts. Resolve them as they arise and continue the rebase process. After resolving conflicts, you can continue the rebase with `git rebase --continue`. Repeat this until the rebase is complete.

6. Push Your Changes
   If you've rebased, you'll need to force-push your changes back to your repository (since you've rewritten history):
   ```bash
   git push origin main --force
   ```

   If you merged, a regular push will suffice:
   ```bash
   git push origin main
   ```

- After completing these steps, your local branch should be synchronized with the upstream `main` branch. However, if you have local commits that are not in the upstream, they will still be present in your local branch.
- If there are conflicts during the `git merge` step, you will need to resolve these conflicts manually before you can complete the merge.
- It's a good practice to regularly fetch and merge changes from the upstream, especially before starting new work or pushing changes to your fork. This minimizes conflicts and keeps your fork up-to-date.

By following these steps, you can easily keep track of any updates or changes that have been made to the central repository and ensure that your fork stays up-to-date with those changes.

### Pull Request (PR)
This is something beyond the scope of this lab. However, I would like to introduce this concept to you. You can read more about PR [here](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests).

Here's a breakdown of how pull requests work and their impact:

1. **Scope of a Pull Request:**
   - A pull request is essentially a request to merge one branch into another within the same repository or across forks.
   - When you create a pull request, it is isolated to your own fork and the branch from which you're making the request.

2. **Impact on Your Repository:**
   - Within your repository, the pull request tracks the changes between the base branch (the one you want to merge into) and the compare branch (the one with your changes).
   - The PR allows for discussion, review, and additional commits until it gets merged or closed.

3. **No Direct Effect on Others' Repositories:**
   - Your pull request does not directly affect the state of anyone else's repository until it is merged by a maintainer of the target repository.
   - Other contributors will not see your changes in their own forks or clones unless they explicitly fetch and merge your branch into their own branches.

4. **After Merging:**
   - If your pull request is accepted and merged into the main repository (the upstream repository from which you forked), then your changes become part of that main repository.
   - Once merged, other contributors can pull the changes from the main repository into their own forks or clones.

5. **Best Practices:**
   - Keep your PRs focused and small to facilitate easier review and quicker integration.
   - Regularly sync your fork with the main repository to avoid merge conflicts.
   - Engage in any discussion or requests for changes that arise during the PR review process.

6. **Open PRs and Syncing:**
   - If you're syncing your fork with the upstream repository while having an open pull request, ensure that your syncing process (merge or rebase) does not disrupt the PR. Typically, changes to other branches or unrelated commits won't affect the PR.

In summary, your pull requests are confined to your own repository and the upstream repository you wish to merge your changes into. They provide a controlled environment for proposing changes, receiving feedback, and collaboratively refining your contributions.
