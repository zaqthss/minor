import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import rand_score

repeat_time = 3

def mainKmeans(dataname, method):
    totalRI = 0.0
    base_path_clean = '../data/clean/Cluster/'
    train_file_clean = base_path_clean + dataname + '/' + dataname + '_TRAIN.csv'
    # the first rows are labels
    df_train = pd.read_csv(train_file_clean, header=None)
    train_y = df_train.iloc[0, 1:]
    n_class = train_y.nunique()
    if method == 'clean':
        df_train = pd.read_csv(train_file_clean, header=None)
        train_x = df_train.iloc[1:, 1:].transpose()
        kmeans = KMeans(n_clusters=n_class, random_state=3)
        kmeans.fit(train_x)
        labels = kmeans.labels_
        totalRI += rand_score(train_y, labels)
    else:
        if model == 'Dirty':
            base_path = '../data/dirty/Cluster/'
        else:
            base_path = '../data/repaired/Cluster/'
        # repaired data do not have labels in the first row
        for i in range(repeat_time):
            seed = str(i)
            train_file = base_path + dataname + '/' + seed + '_' + dataname + '_' + method + '_TRAIN.csv'
            df_train = pd.read_csv(train_file, header=None)
            train_x = df_train.iloc[:, 1:].transpose()
            kmeans = KMeans(n_clusters=n_class, random_state=3)
            kmeans.fit(train_x)
            labels = kmeans.labels_
            totalRI += rand_score(train_y, labels)
        totalRI /= repeat_time
    return totalRI


if __name__ == '__main__':
    data_names = ['DistalPhalanxTW', 'EOGVerticalSignal']
    models = ['clean', 'Dirty', 'MINOR-B', 'MINOR-U', 'IMR', 'MTCSC', 'VARX', 'Akane']
    for data_name in data_names:
        print(data_name + ' Cluster RI: ')
        for model in models:
            result = mainKmeans(data_name, model)
            print(model + ': ' + str(result))
        print('\n')
