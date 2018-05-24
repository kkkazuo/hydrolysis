import random
from datetime import datetime

_source_str = 'abcdefghijklmnopqrstuvwxyz0123456789'


def issue_id():
    random_num = ''.join([random.choice(_source_str) for _ in range(7)])
    today = datetime.today()
    date_str = str(datetime.date(today))
    return date_str+random_num


if __name__ == "__main__":
    print(issue_id())
