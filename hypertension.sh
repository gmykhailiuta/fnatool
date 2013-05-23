grep [hH]ypertension s*.hea | grep -v No | cut -d: -f1 | sort | uniq | cut -d. -f1
