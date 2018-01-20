# Version 0.1

## Round 1

### Test environments
* local OS X install, R 3.4.1
* ubuntu 12.04 (on travis-ci), R 3.4.1
* win-builder (devel and release)

### R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

### Reverse dependencies

This is a new release, so there are no reverse dependencies.

### Reviewer Comments

Swetlana Herbrant 09-20-17:

Thanks, can you provide a reference in the 'Description' field of your DESCRIPTION file in the form authors (year) <doi:...> or <arXiv:...> (no space after 'doi:' and 'arXiv:')?

Any reason why 'Escalation With Overdose Control' is capitalized?

Your examples are wrapped in \dontrun{}, hence nothing gets tested. 
Please unwrap the examples if that is feasible and if they can be executed in < 5 sec for each Rd file or create additionally small toy examples. Something like \examples{
        examples for users:
        executable in < 5 sec
        for checks
        \dontshow{
               examples for checks:
               executable in < 5 sec together with the examples above
               not shown to users
        }
        donttest{
               further examples for users (not used for checks)
        }
        \dontrun{
               examples with e.g. very long run times
               (not used for checks)
        }
}
would be desirable.


## Round 2

### Submission comments

Addressed all previous comments.

The correct DOI in DESCRIPTION is 10.1002/(SICI)1097-0258(19980530)17:10<1103::AID-SIM793>3.0.CO;2-9, but R Check considers only 10.1002/(SICI)1097-0258(19980530)17:10<1103::AID-SIM793 ignoring the
last part >3.0.CO;2-9 .

### Reviewer Comments

Swetlana Herbrant 10-04-17:

Thanks, we get:
Found the following (possibly) invalid DOIs:
   DOI: 10.1002/(SICI)1097-0258(19980530)17:10<1103::AID-SIM793
     From: DESCRIPTION
     Status: Not Found
     Message: 404

Please do not capitalize 'Escalation with Overdose Control' in your description.

Please add space after periods.

Please add more small executable examples in your Rd-files for checks and users.


Please fix and resubmit.

Uwe Ligges 10-04-17:

Reason: the DOI contains the special chars "<" and ">", hence please write it as


URLencode("10.1002/(SICI)1097-0258(19980530)17:10<1103::AID-SIM793>3.0.CO;2-9")

i.e.<doi:10.1002/(SICI)1097-0258(19980530)17:10%3C1103::AID-SIM793%3E3.0.CO;2-9>.

## Round 3

### Submission comments

'Escalation with Overdose Control' is no longer capitalized and space after 
periods was added.

DOI written in the form  URLencode("10.1002/(SICI)1097-0258(19980530)17:10<1103::AID-SIM793>3.0.CO;2-9")

Simple examples were added to all Rd-files.

### Reviewer comments

Thanks, please write the DOI in the form:

<doi:10.1002/(SICI)1097-0258(19980530)17:10%3C1103::AID-SIM793%3E3.0.CO;2-9>

Please fix and resubmit.

### Round 4

### Submission comments

DOI written in the form 
<doi:10.1002/(SICI)1097-0258(19980530)17:10%3C1103::AID-SIM793%3E3.0.CO;2-9>

---

# Version 0.1.1

### Test environments
* local OS X install, R 3.4.1
* ubuntu 12.04 (on travis-ci), R 3.4.1
* win-builder (devel and release)

### R CMD check results

R CMD check results
0 errors | 0 warnings | 0 notes

### Reverse dependencies

This is a new release, so there are no reverse dependencies.
