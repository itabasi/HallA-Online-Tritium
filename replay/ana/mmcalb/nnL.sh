#! /bin/sh

days="0916"

H1="./bin/ana_Lambda -f ../run_list/nnlambda/Lambda_small_H1.list -p param/f1_tuned_Lambda_twc.param -H1 -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/Lambda_small_H1_$days.root -w ../pdf/mmass/ana_Lambda/Lambda_small_H1_$days.pdf"

H2="./bin/ana_Lambda -f ../run_list/nnlambda/Lambda_small_H2.list -p param/f1_Lambda_phase2_tuned.param -H1 -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/Lambda_small_H2_$days.root -w ../pdf/mmass/ana_Lambda/Lambda_small_H2_$days.pdf"


hadd_H="hadd ../rootfiles/mmass/ana_Lambda/Lambda_small_H_$days.root ../rootfiles/mmass/ana_Lambda/Lambda_small_H1_$days.root ../rootfiles/mmass/ana_Lambda/Lambda_small_H2_$days.root"

H3="./bin/ana_Lambda -f ../run_list/nnlambda/Lambda_small_T.list -p param/f1_Lambda_phase2_tuned.param -H2 -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/Lambda_small_T_$days.root -w ../pdf/mmass/ana_Lambda/Lambda_small_T_$days.pdf"

T1="./bin/ana_Lambda -f ../run_list/nnlambda/nnL_small_1.list -p param/f1_tuned_Lambda_twc.param -T -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/nnL_small_1_$days.root -w ../pdf/mmass/ana_Lambda/nnL_small_1_$days.pdf"

T2="./bin/ana_Lambda -f ../run_list/nnlambda/nnL_small_234.list -p param/f1_Lambda_phase2_tuned.param -T -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/nnL_small_2_$days.root -w ../pdf/mmass/ana_Lambda/nnL_small_2_$days.pdf"

hadd_T="hadd ../rootfiles/mmass/ana_Lambda/nnL_small_$days.root ../rootfiles/mmass/ana_Lambda/nnL_small_1_$days.root ../rootfiles/mmass/ana_Lambda/nnL_small_2_$days.root"


eval $H1
eval $H2
eval $H3
eval $T1
eval $T2
eval $hadd_H
eval $hadd_T
