import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import f1_score

repeat_time = 3


def mainKnn(K, dataname, method):
    totalF1 = 0.0
    base_path_clean = '../data/clean/Classification/';
    train_file_clean = base_path_clean + dataname + '/' + dataname + '_TRAIN.csv';
    test_file_clean = base_path_clean + dataname + '/' + dataname + '_TEST.csv';
    # the first rows are labels
    df_train = pd.read_csv(train_file_clean, header=None)
    train_y = df_train.iloc[0, 1:]
    df_test = pd.read_csv(test_file_clean, header=None)
    test_y = df_test.iloc[0, 1:]
    if method == 'clean':
        df_train = pd.read_csv(train_file_clean, header=None)
        train_x = df_train.iloc[1:, 1:].transpose()
        df_test = pd.read_csv(test_file_clean, header=None)
        test_x = df_test.iloc[1:, 1:].transpose()
        knn = KNeighborsClassifier(n_neighbors=K)
        knn.fit(train_x, train_y)
        prediction = knn.predict(test_x)
        totalF1 = f1_score(test_y, prediction, average='macro')
    else:
        if model == 'Dirty':
            base_path = '../data/dirty/Classification/'
        else:
            base_path = '../data/repaired/Classification/'
        # repaired data do not have labels in the first row
        for i in range(repeat_time):
            seed = str(i)
            train_file = base_path + dataname + '/' + seed + '_' + dataname + '_' + method + '_TRAIN.csv'
            test_file = base_path + dataname + '/' + seed + '_' + dataname + '_' + method + '_TEST.csv'
            df_train = pd.read_csv(train_file, header=None)
            train_x = df_train.iloc[:, 1:].transpose()
            df_test = pd.read_csv(test_file, header=None)
            test_x = df_test.iloc[:, 1:].transpose()
            knn = KNeighborsClassifier(n_neighbors=K)
            knn.fit(train_x, train_y)
            prediction = knn.predict(test_x)
            totalF1 += f1_score(test_y, prediction, average='macro')
        totalF1 /= repeat_time
    return totalF1


if __name__ == '__main__':
    data_names = ['Car', 'Lightning2']
    Ks = [5, 5]
    models = ['clean', 'Dirty', 'MINOR-B', 'MINOR-U', 'IMR', 'MTCSC', 'VARX', 'Akane']
    for i in range(len(data_names)):
        print(data_names[i] + ' Classification F1-scores: ')
        for model in models:
            result = mainKnn(Ks[i], data_names[i], model)
            print(model + ': ' + str(result))
        print('\n')
